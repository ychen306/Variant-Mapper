from __future__ import division
from flask import Flask, Response, request, render_template, session
import ijson
import requests
from gevent.pool import Pool
from gevent import monkey
import simplejson as json
import collections
import re
from functools import partial
from settings import GOOGLE_API_KEY
from models import CachedReads
from variant_mapper import match, get_ref_length

monkey.patch_all()

app = Flask(__name__)
app.secret_key = '3B69607696501244C5516219551601DE0700FAD017B62E75527C2EA70C740C35'

REPOSITORIES = {
    'google': 'https://www.googleapis.com/genomics/v1beta/',
    'ncbi': 'http://trace.ncbi.nlm.nih.gov/Traces/gg/',
    'ebi': 'http://193.62.52.16/'
}
MIN_DISTANCE = 100
Coordinate = collections.namedtuple('Coordinate', 'chrom start end')


class HTTPStream():
    '''
    Turn requests.Response into a stream
    '''
    def __init__(self, response):
        self.response = response
    def read(self, size):
        if not hasattr(self, 'iter'):
            self.iter = iter(self.response.iter_content(chunk_size=size))
        try:
            return next(self.iter)
        except StopIteration:
            return ''



def find_cached_read(report, repository, readset_id):
    '''
    Find cached read by report 
    '''
    return CachedReads.find_one({
            'repository': repository,
            'readsetId': readset_id,
            'referenceSequenceName': report['chrom'],
            'start': {'$lte': report['seqStart']},
            'end': {'$gt': report['seqEnd']}
        })


def push_coordinates(coordinates, new_coord):
    '''
    Given a list of existing coordinates and a single new coordinate, return a new list of coordinates,
    optimizing search operation of reads.
    '''
    for i, coord in enumerate(coordinates):
        if coord.chrom != new_coord.chrom:
            continue
        #calculate size of the distance between two coordinates
        upper = max(coord.start, new_coord.start)
        lower = min(coord.end, new_coord.end)
        coord_index = i
        #meaning the two coordinates overlap
        if upper <= lower or (upper - lower) <= MIN_DISTANCE:
            start = min(coord.start, new_coord.start)
            end = max(coord.end, new_coord.end)
            coordinates[i] = Coordinate(coord.chrom, start, end)
            break
    else:
        coordinates.append(new_coord)
        coord_index = len(coordinates) - 1
    #coord_index tells the index of the newly inserted coordinate
    return coordinates, coord_index


def make_read_search(repo_id, readset_id, coordinate, coord_index):
    '''
    Construct a search on read api, return the a requests.request object
    '''
    search_data = {
        'readsetIds': [readset_id],
        'sequenceName': coordinate.chrom,
        'sequenceStart': coordinate.start,
        'sequenceEnd': coordinate.end
    }
    if repo_id == 'google':
        api_key = 'key=%s&' % GOOGLE_API_KEY
    else:
        api_key = ''
    search_url = '%sreads/search?%scoord_index=%s' % (REPOSITORIES[repo_id], api_key, coord_index)
    return requests.post(search_url,
                        data=json.dumps(search_data),
                        stream=True,
                        headers={'Content-Type': 'application/json; charset=UTF-8'})




@app.route('/<repo_id>/<path:endpoint>', methods=['GET', 'POST'])
def ga_api(repo_id, endpoint):
    '''
    Makes GA4GH read api call for the front end.
    '''
    url = REPOSITORIES[repo_id] + endpoint
    if repo_id == 'google':
        url += '?key=%s' % GOOGLE_API_KEY
    url += '&coord_id=1'
    options = {'stream': True}
    if request.method == 'POST':
        options['data'] = request.data
        options['headers'] = {'Content-Type': 'application/json; charset=UTF-8'}
    response = requests.request(request.method, url, **options)
    return Response(response.iter_content(),
                        content_type='application/json; charset=UTF-8',
                        status=response.status_code)



@app.route('/match_reports', methods=['POST'])
def match_reports():
    report_set = json.loads(request.data)
    coordinates = []
    coord_indices = {}
    matched_reports = []

    for report in report_set['reports']:
        if report['clinicalSignificance'] in ('Uncertain significance', 'not provided', 'conflicting data from submitters', 'other'):
            continue
        new_coord = Coordinate(report['chrom'], report['seqStart'], report['seqEnd'])
        coordinates, coord_index = push_coordinates(coordinates, new_coord)
        cached_read = find_cached_read(report, report_set['repository'], report_set['readsetId'])
        if cached_read:
            if match(report, cached_read):
                matched_reports.append(report)
        else:
            #look it up with read search API
            coord_indices.setdefault(coord_index, []).append(report)

    read_search = partial(make_read_search,
                        report_set['repository'],
                        report_set['readsetId'])

    def search_and_match(coord_index):
        coord = coordinates[coord_index]
        result = read_search(coord, coord_index)
        reads = ijson.items(HTTPStream(result), 'reads.item')
        reports = coord_indices[coord_index]
        reports_visited = []
        for read in reads:
            read_coord = Coordinate(read['referenceSequenceName'],
                                    read['position'],
                                    read['position']+get_ref_length(read['cigar']))
            covered_reports = [report for report in reports
                                    if (read_coord.chrom == report['chrom'] and
                                        read_coord.start <= report['seqStart'] and
                                        read_coord.end > report['seqEnd'])]
            if covered_reports:
                new_matched_reports = [report for report in covered_reports 
                                        if (report['reportId'] not in reports_visited and
                                            match(report, read))]
                matched_reports.extend(new_matched_reports)
                reports_visited.extend([report['reportId'] for report in covered_reports
                                            if report['reportId'] not in reports_visited])
                #push read into cache:
                read['repository'] = report_set['repository']
                read['readsetId'] = report_set['readsetId']
                read['start'] = read['position']
                read['end'] = read['position'] + get_ref_length(read['cigar'])
                CachedReads.save(read)
            if len(reports_visited) >= len(reports):
                break

    pool = Pool(50)
    pool.imap_unordered(search_and_match, coord_indices)
    pool.join()

    return Response(json.dumps(matched_reports), content_type='application/json; charset=UTF-8')




@app.route('/register')
def register():
    session['repo'] = request.args['repo']
    session['dataset'] = request.args['dataset']
    session['readset'] = request.args['readset']
    return Response()



@app.route('/', methods=['GET', 'POST'])
def index():
    app_init = ''
    if session.get('repo'):
        app_init += 'currentRepo="%s"; '% session['repo']
        if session.get('dataset'):
            app_init += 'currentDataset="%s"; '% session['dataset']
        if session.get('readset'):
            app_init += 'currentReadset=%s; '% session['readset']

    if request.form.get('term'):
        searchTerms = 'mapTerm("%s");'% request.form['term'];
    else:
        searchTerms = '';

    return render_template('index.html', app_init=app_init, searchTerms=searchTerms, term=request.form.get('term', ''))




if __name__ == '__main__':
    app.run(debug=True, port=7000)