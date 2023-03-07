from typing import Dict


def check_version(json_cats: Dict, version: str) -> None:
    """
    Check the version.
    Sends an exception if the version is different from the expected version.
    """
    schema_version = json_cats['metaData']['schemaVersion']
    if schema_version != version:
        raise ValueError(f"invalid version: {schema_version}")


def iter_variants(json_cats: Dict, key: str):
    """Extract Mutation"""
    variants_set: Dict = json_cats.get('variants', {})
    for variant in variants_set.get(key, []):
        yield variant
