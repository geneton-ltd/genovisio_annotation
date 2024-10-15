import argparse
import gzip
import json

from annotation.src import cnv_region, collections, core, schemas

import enum

class _CollectionName(enum.StrEnum):
    BENIGN_CNV = 'Benign_CNV'
    BENIGN_CNV_GS_INNER = 'Benign_CNV_GS_inner'
    BENIGN_CNV_GS_OUTER = 'Benign_CNV_GS_outer'
    REGULATORY = 'Regulatory'
    GNOMAD = 'GnomAD'
    HI_GENE = 'HI_gene'
    HI_REGION = 'HI_region'
    GENES = 'Genes'

    def values(self):
        return [item.value for item in self.__class__]


def get_whole_annotation(input_str: str, mongodb_uri: str, db_name: str) -> schemas.Annotation:
    region = cnv_region.CNVRegion.build_from_string(input_str)
    genovisio_db = collections.GenovisioDB(mongodb_uri, db_name)
    results = {collection_name: genovisio_db.find_intersections(collection_name, region) for collection_name in _CollectionName}
    cnv_annotation = schemas.CNVRegionAnnotation(
        chr=region.chr,
        start=region.start,
        end=region.end,
        cnv_type=region.cnv_type,
        length=region.length,
        name=region.name,
        cytogenetic_position=region.cytogenetic_position
    )

    return schemas.Annotation(
        cnv=cnv_annotation,
        benign_cnv=results[_CollectionName.BENIGN_CNV],
        benign_cnv_gs_inner=results[_CollectionName.BENIGN_CNV_GS_INNER],
        benign_cnv_gs_outer=results[_CollectionName.BENIGN_CNV_GS_OUTER],
        regulatory=results[_CollectionName.REGULATORY],
        gnomad=results[_CollectionName.GNOMAD],
        hi_gene=results[_CollectionName.HI_GENE],
        hi_region=results[_CollectionName.HI_REGION],
        genes=results[_CollectionName.GENES]
    )


def main():
    parser = argparse.ArgumentParser(description='Classify CNV and/or find intersecting items in MongoDB collections.')
    parser.add_argument('input', help='Input string in the form "chr1:10000-20000/del". CNV type should be del/dup/loss/gain.')
    parser.add_argument('--output', '-o', help='Filename, where to store the resulting json. (add .gz for gzipped output)', required=True)
    parser.add_argument('--mongodb_uri', help='MongoDB full URI', default='mongodb://localhost:27017/')
    parser.add_argument('--db_name', help='MongoDB database name', default='genovisio')
    args = parser.parse_args()


    annotation = get_whole_annotation(
        input_str=args.input,
        mongodb_uri=args.mongodb_uri,
        db_name=args.db_name
    )

    if args.output.endswith('.gz'):
        with gzip.open(args.output, 'wt', encoding=core.ENCODING) as f:
            f.write(json.dumps(annotation, default=str, indent=core.JSON_INDENT_LEVEL))
    else:
        with open(args.output, 'w', encoding=core.ENCODING) as f:
            f.write(json.dumps(annotation, default=str, indent=core.JSON_INDENT_LEVEL))


if __name__ == '__main__':
    main()
