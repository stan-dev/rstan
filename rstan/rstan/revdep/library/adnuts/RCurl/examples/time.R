library(RCurl)

body = '<?xml version="1.0" encoding="utf-8"?>
<soap:Envelope xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:soap="http://schemas.xmlsoap.org/soap/envelope/">
  <soap:Body>
    <getTimeZoneTime xmlns="http://www.Nanonull.com/TimeService/">
      <timezone>GMT</timezone>
    </getTimeZoneTime>
  </soap:Body>
</soap:Envelope>'

curlPerform( url = "http://www.nanonull.com/TimeService/TimeService.asmx",
             httpheader = c(accept = "multipart/*",
                            "Content-Type" = "text/xml; charset=utf-8",
                            SOAPAction = '"http://www.Nanonull.com/TimeService/getTimeZoneTime"'),
             postfields = body)
