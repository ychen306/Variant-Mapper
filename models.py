from pymongo import MongoClient
from settings import MONGO_URL, MONGO_DB

db = MongoClient(MONGO_URL)[MONGO_DB]
CachedReads = db['read-cache']
'''
Schema for CachedReads
{
    _id: ObjectId,
    repository: string,
    start: int,
    end: int,
    ...GA4GH read...

}
'''
