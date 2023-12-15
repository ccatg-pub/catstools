{
    "metaData": {
        "schemaVersion": "1.3.0",
        "referenceGenome": {
            "name": "GRCh37",
            "grcRelease": "GRCh37",
            "descriptions": [
                "sample descriptions"
            ]
        },
        "comments": []
    },
    "testInfo": {
        "testId": "20201122112233",
        "testType": "tumor-only (cell-free)",
        "softwareName": "softwareName",
        "softwareVersion": "softwareVersion",
        "panelName": "TestSamplePanel",
        "panelVersion": "XYZ"
    },
    "otherBiomarkers": [
        {
            "itemId": "biomarker-1",
            "biomarkerType": "MSI",
            "state": "null",
            "sampleItemId": "sequence-1-tumor-dna",
            "reported": true
        },
        {
            "itemId": "biomarker-2",
            "biomarkerType": "TMB",
            "biomarkerMetrics": [
                {
                    "value": 1.21,
                    "unit": "mutations-per-megabase",
                    "type": "mutations-per-megabase"
                }
            ],
            "state": "low",
            "sampleItemId": "sequence-1-tumor-dna",
            "reported": true
        }
    ],
    "sequencingSamples": [
        {
            "itemId": "sequence-1-tumor-dna",
            "tumorOrNormal": "tumor",
            "nucleicAcid": "DNA"
        }
    ]
}