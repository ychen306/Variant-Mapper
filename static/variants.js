datasets = {
    'Google': {
        '1000 Genomes': '10473108253681171589',
        'DREAM SMC Challenge': '337315832689',
        'PGP': '383928317087',
        'Simons Foundation' : '461916304629'
    },
    'NCBI': {
        'SRP034507': 'SRP034507',
        'SRP029392': 'SRP029392'
    },
    'EBI': {
        'All data': 'data'
    }
};



var variantMapper = function(repo, readsetId, $scope) {

    var totalReportSets = 0,
        processedReportSets = 0,
        matchesSubmitted = 0;
    $.getXML = function(url, callback) {
        $.ajax({
            type: 'GET',
            url: url,
            dataType: 'xml',
            success: callback
        });
    };

    var plus = function(a, b) {
        return a + b;
    };

    var minus = function(a, b) {
        return a - b;
    };


    var evalArith = function(equation) {
        /*
        Evaluate plus and minus equation to deal with weird HGVS notation
        */
        var result = undefined;
        if (equation.indexOf('+') !== -1) {
            var operation = plus,
                operants = equation.split('+');
        } else if (operation.index('-') !== -1) {
            var operation = minus,
                operants = equation.split('-');
        } else {
            result = parseInt(equation);
        }
        if (!result) {
            result = operation(parseInt(operants[0]), parseInt(operants[1]));
        }
        return result;
    };


    var displayReports = function(reports) {
        /*
        Render matched reports.
        */
        processedReportSets += 1;
        console.log(processedReportSets + '/' + totalReportSets);
        reports.forEach(function(report){
            $scope.$apply(function() {
                $scope.matchedReports.push(report);
            });
        });
        $scope.$apply(function() {
            if ($scope.stage === "Matching reports with readset") {
                $scope.progress = (processedReportSets/totalReportSets*100).toFixed(2) + "%";
            }
            if (processedReportSets >= totalReportSets) {
                $scope.loadingReports = false;
            }
        });
    };

    var DocumentParsingError = function() {
        this.message = 'Error parsing document from clinvar';
        this.name = 'DocumentParsingError';
    };


    var constructClinvarReport = function(doc){
        var hgvsPattern = /.*?[cg]\.([\d_\+\-\?]+)(\S+)(?: \(.+?\))?/,
            report = new Object(),
            hgvs = doc.find('variation_set_name').first().text().match(hgvsPattern);
        var seqRange = hgvs[1];
        report.chrom = doc.find('chr_sort').first().text();
        report.chrom = report.chrom.substring(report.chrom.lastIndexOf('0') + 1);
        report.seqStart = parseInt(doc.find('location_sort').first().text());
        report.variation = hgvs[2];
        report.wholeVariation = hgvs[0];
        report.strand = doc.find('strand').first().text();
        if (seqRange.indexOf('_') === -1) {
            report.seqEnd = report.seqStart;
        } else {
            var seqRange = seqRange.split('_');
            var start = evalArith(seqRange[0]),
                end = evalArith(seqRange[1]);
            report.seqEnd = report.seqStart + end - start;
        }
        report.gene = doc.find('gene_sort').first().text();
        report.evaluated = doc.find('last_evaluated').first().text();
        report.trait = doc.find('trait_name').first().text();
        report.clinicalSignificance = doc.find('description').first().text();
        report.reviewStatus = doc.find('review_status').first().text();
        report.reportId = doc.attr('uid');
        report.title = doc.find('title').first().text();
        for (var key in report) {
            if (!report[key]) {
                throw DocumentParsingError();
            }
        }
        return report
    }



    var submitMatch = function(reportSet) {
        /*
        Submit a set of reports to server, which will then check if any of the reports matches the read set.
        */
        $.ajax({
            type: 'POST',
            url: '/match_reports',
            beforeSend: function(request) {
                request.setRequestHeader('Content-Type', 'application/json; charset=UTF-8');
            },
            data: JSON.stringify(reportSet),
            success: displayReports,
            error: function() {
                setTimeout(function() {
                    submitMatch(reportSet);
                }, 10000);
            }
        });
    };


    var getClinvarReports = function(idList) {
        var summaryUrl = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=';
        $(idList).find('Id').each(function() {
            summaryUrl += ($(this).text() + ',');
        });
        $.getXML(summaryUrl, function(xml) {
            var reportSet = new Object();
            reportSet.reports = new Array();
            reportSet.repository = repo.toLowerCase();
            reportSet.readsetId = readsetId;
            $(xml).find('DocumentSummary').each(function(){
                try {
                    var report = constructClinvarReport($(this));
                } catch (e) {
                    return;
                }
                reportSet.reports.push(report);
            });
            if (reportSet.reports.length) {
                submitMatch(reportSet);
            } else {
                processedReportSets += 1;
            }
            matchesSubmitted += 1;
        });
    };




    return {
        map: function(term) {
            $scope.matchedReports = new Array();
            var retMax = 50,
                searchUrl = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&retmax='+retMax+'&term='+term,
                trialUrl = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&retmax=0&term='+term;
            $.getXML(trialUrl, function(initSearch) {
                var count = parseInt($(initSearch).find('Count').first().text()),
                    retStart = 0;
                if (count > 0) {
                    function searchReports() {
                        totalReportSets += 1;
                        $.getXML(searchUrl+'&retstart='+retStart, getClinvarReports);
                        retStart += retMax;
                        if (retStart < count) {
                            setTimeout(searchReports, 5);
                            $scope.$apply(function() {
                                $scope.progress = (retStart/count*100).toFixed(2) + "%";
                            });
                        } else {
                            $scope.$apply(function() {
                                $scope.progress = "100%";
                                setTimeout(function() {
                                    $scope.progress = "0%";
                                    $scope.stage = "Matching reports with readset";
                                }, 800);
                            });
                        }
                    }
                    searchReports();
                } else {
                    //not reports found in clinvar.
                    $scope.$apply(function() {
                        $scope.loadingReports = false;
                    });
                }
            });
        }
    };
};


var getReadsets = function(repository, datasetName, $scope, pager) {
    var dataset = datasets[repository][datasetName],
        searchData = {
            datasetIds: [dataset]
        },
        requestData = {
            method: 'POST',
            url: repository.toLowerCase()+'/'+'readsets/search',
            headers: {'Content-Type': 'application/json; UTF-8'}
        };
    if (typeof pager !== 'object') {
        pager = {};
    }
    if (pager.last === true) {
        var tokenIndex = $scope.pageTokens.length - 3;
        searchData.pageToken = $scope.pageTokens[tokenIndex];
        $scope.pageTokens.pop();
        $scope.pageTokens.pop();
    } else if (pager.next === true) {
        var tokenIndex = $scope.pageTokens.length - 1;
        searchData.pageToken = $scope.pageTokens[tokenIndex];
    }
    $scope.pageTokens.push(undefined); //placeholder.
    requestData.body = JSON.stringify(searchData);
    //remove this after NCBI changes pageToken to nextPageToken
    if (repository == 'NCBI') {
        var nextPageTokenName = 'pageToken';
    } else {
        var nextPageTokenName = 'nextPageToken';
    }
    $scope.loading = true;
    oboe(requestData)
        .node('readsets.*', function(readset){
            if ($scope.currentDataset != datasetName || $scope.currentRepo != repository) {
                this.forget();
                return;
            } 
            try {
                species = readset.fileData[0].refSequences[0].species;
                if (!species) {
                    species = "Unknown";
                }
                readset.species = species;
            } catch (e) {}
            $scope.$apply(function() {
                $scope.readsets.push(readset);
            });
        })
        .node(nextPageTokenName, function(nextPageToken) {
            $scope.$apply(function() {
                $scope.pageTokens[$scope.pageTokens.length - 1] = nextPageToken; //replace undefined with real token.
            });
        });

};
