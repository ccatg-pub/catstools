from typing import Dict

from catstools.formatupdate.utils import check_version


def update_meta(json_cats: Dict):
    """Update metaData"""
    json_cats['metaData']['schemaVersion'] = "1.0.1"


def run(json_cats: Dict):
    check_version(json_cats, "1.0")
    update_meta(json_cats)
