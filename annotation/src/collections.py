from dataclasses import dataclass, field
from typing import Any

import pymongo

from annotation.src import cnv_region


@dataclass
class GenovisioDB:
    uri: str
    name: str
    client: pymongo.MongoClient = field(init=False)
    db: pymongo.database.Database = field(init=False)

    def __post_init__(self) -> None:
        self.client = pymongo.MongoClient(self.uri)
        self.db = self.client[self.name]

    def find_intersections(self, collection_name: str, region: cnv_region.CNVRegion) -> list[Any]:
        """
        Find intersecting items within a single collection using indexed queries.
        """
        # Use indexed query to find potential intersections
        query = {"chromosome": region.chr, "start": {"$lte": region.end}, "end": {"$gte": region.start}}
        return list(self.db[collection_name].find(query))
