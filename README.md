Variant Mapper specifically written for the Read API of Global Alliance for Genomics and Health. This web application searches for reports filed by physicians and researches and compare those reports with a readset exposed by the Read API. A live server is running [here](http://192.241.244.189:7000). 
## Running the server

### Prerequisites
* MongoDB
* Install required Python packages with following command:
```
$ pip install -r requirements.txt
```

### Configuration
* Edit `settings.py.defalut` and set your own Mongodb settings.
* Set your own `SECRET_KEY`.
* You will also need a Google API Key.
* Save `settings.py.default` as `settings.py`

### Run the server
* To run the server locally, use command below:
```
$ python main.py
```