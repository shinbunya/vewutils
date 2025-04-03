def getQH(siteids,staDT,endDT,cachedir=None):
    import requests
    import urllib.parse
    from datetime import datetime
    import time
    import hashlib
    import os
    from pathlib import Path
    import json

    if type(siteids) != list:
        siteids = [siteids]

    stStaDT = staDT.strftime("%Y-%m-%dT%H:%M%z")
    stEndDT = endDT.strftime("%Y-%m-%dT%H:%M%z")

    timqs = []
    valqs = []
    timhs = []
    valhs = []

    for siteid in siteids:
        params = {'format': 'json', 'sites': siteid, 'startDT': stStaDT, 'endDT': stEndDT, 'parameterCd': '00060,00065', 'siteType': 'ST', 'siteStatus':'all'}
        params_encoded = urllib.parse.urlencode(params)
        api_url = 'https://nwis.waterservices.usgs.gov/nwis/iv/?'+params_encoded
        if cachedir:
            hasher = hashlib.md5()
            hasher.update(api_url.encode('utf-8'))
            cachefile = os.path.join(cachedir, hasher.hexdigest()+'.json')
        else:
            cachefile = None
        print('Target = {}'.format(api_url))
        if cachefile and Path(cachefile).is_file():
            print('Found in local cache: {:s}. Loading.'.format(cachefile))
            with open(cachefile, 'r') as f:
                json_res = json.load(f)
        else:
            print('Not found in local cache. Downloading.')
            cnt = 0
            while True:
                try:
                    response = requests.get(api_url)
                    json_res = response.json()
                    break
                except requests.JSONDecodeError as e:
                    cnt += 1
                    print('JSONDecodeError: {}'.format(e))
                    print('Response content: {}'.format(response.content))
                    if cnt < 10:
                        print('Retrying in 2 seconds...')
                    else:
                        print('Failed to decode JSON response after 10 attempts. Raising exception.')
                        raise e
                time.sleep(2)
            


            if cachefile:
                with open(cachefile, 'w') as f:
                    json.dump(json_res, f)

        # Discharge
        nQ = len(json_res['value']['timeSeries'][0]['values'][0]['value'])
        timq = [datetime.strptime(json_res['value']['timeSeries'][0]['values'][0]['value'][i]['dateTime'],'%Y-%m-%dT%H:%M:%S.%f%z') for i in range(nQ)]
        valq = [float(json_res['value']['timeSeries'][0]['values'][0]['value'][i]['value'])*0.3048**3 for i in range(nQ)]

        # Gage height
        nH = len(json_res['value']['timeSeries'][1]['values'][0]['value'])
        timh = [datetime.strptime(json_res['value']['timeSeries'][1]['values'][0]['value'][i]['dateTime'],'%Y-%m-%dT%H:%M:%S.%f%z') for i in range(nH)]
        valh = [float(json_res['value']['timeSeries'][1]['values'][0]['value'][i]['value'])*0.3048 for i in range(nH)]

        # Store
        timqs.append(timq)
        valqs.append(valq)
        timhs.append(timh)
        valhs.append(valh)

    return timqs, valqs, timhs, valhs