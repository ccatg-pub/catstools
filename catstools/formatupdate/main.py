from collections import OrderedDict
from itertools import islice
from typing import Callable, Dict, Iterator, Optional

from catstools.formatupdate import v1_0_0, v1_0_1, v1_1_0, v1_2_0, v1_2_1

VERSION_MAP = OrderedDict([
    ("1.0", v1_0_0.run),
    ("1.0.1", v1_0_1.run),
    ("1.1.0", v1_1_0.run),
    ("1.2.0", v1_2_0.run),
    ("1.2.1", v1_2_1.run),
])
VERSION_LIST = ", ".join(list(VERSION_MAP.keys()))


def get_version(json_cats: Dict) -> str:
    """Get cats format version."""
    return json_cats['metaData']['schemaVersion']


def iter_tasks(from_version: str, to_version: Optional[str]) -> Iterator[Callable]:
    """Build the tasks required for the update."""
    versions = list(VERSION_MAP.keys())
    to_index = versions.index(to_version) if to_version else None
    yield from islice(VERSION_MAP.values(), versions.index(from_version), to_index)


def run(json_cats: Dict, to_version: Optional[str]) -> bool:
    """Perform updates in cats format."""
    if to_version is None or to_version in VERSION_MAP.keys():
        from_version = get_version(json_cats)
        for task in iter_tasks(from_version, to_version):
            task(json_cats)
        return True
    else:
        print("usage: catstools format_update [-h] --input INPUT --output OUTPUT [--to-version TO_VERSION]")
        print(f"catstools format_update: error: argument --to-version: invalid choice: {to_version}")
        print(f"(choose from [{VERSION_LIST}])")
        return False
