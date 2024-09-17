def getQH(siteids,staDT,endDT):
    import requests
    import urllib.parse
    from datetime import datetime

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
        print(api_url)
        response = requests.get(api_url)
        json_res = response.json()

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