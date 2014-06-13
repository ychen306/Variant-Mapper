Variant Mapper specifically written for GA4GH Read API

## Running the server

### Prerequisites
* MongoDB
* Install required Python packages with following command:
```
$ pip install -r requirements.txt
```

### Configuration
* Edit `settings.py.defalut` and set your own Mongodb settings. You will also need a Google API Key.
* Save `settings.py.default` as `settings.py`

### Run the server
* To run the server locally, use command below:
```
$ python main.py
```