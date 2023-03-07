from typing import Dict

from catstools.formatupdate.utils import check_version


def run(json_cats: Dict):
    check_version(json_cats, "1.2.1")
