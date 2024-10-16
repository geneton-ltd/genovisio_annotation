import annotation
from annotation.src.schemas import CNVRegionAnnotation


def test_cnv(annotation_cnv: CNVRegionAnnotation):
    assert annotation_cnv.chr == "chr1"
    assert annotation_cnv.start == 10000
    assert annotation_cnv.end == 20000
    assert annotation_cnv.cnv_type == "gain"
    assert annotation_cnv.length == 10000
    assert annotation_cnv.name == "test_name"
    assert annotation_cnv.cytogenetic_position == "test"

    assert annotation_cnv.is_duplication
    assert annotation_cnv.genomic_coord == "chr1:10000-20000"


def test_cnv_has_start_in_region(annotation_cnv: CNVRegionAnnotation):
    assert annotation_cnv.has_start_in_region(9000, 11000)
    assert not annotation_cnv.has_start_in_region(11000, 12000)

    assert annotation_cnv.has_end_in_region(19000, 21000)
    assert not annotation_cnv.has_end_in_region(11000, 12000)


def test_cnv_is_overlapping(annotation_cnv: CNVRegionAnnotation):
    region1 = (9000, 11000)
    assert annotation_cnv.is_overlapping(region1[0], region1[1], annotation.enums.OverlapType.all)
    assert not annotation_cnv.is_overlapping(region1[0], region1[1], annotation.enums.OverlapType.inside_only)
    assert annotation_cnv.is_overlapping(region1[0], region1[1], annotation.enums.OverlapType.partial_both)
    assert not annotation_cnv.is_overlapping(region1[0], region1[1], annotation.enums.OverlapType.partial_end)
    assert annotation_cnv.is_overlapping(region1[0], region1[1], annotation.enums.OverlapType.partial_start)
    assert not annotation_cnv.is_overlapping(region1[0], region1[1], annotation.enums.OverlapType.span_whole_only)

    region2 = (9000, 21000)
    assert annotation_cnv.is_overlapping(region2[0], region2[1], annotation.enums.OverlapType.all)
    assert not annotation_cnv.is_overlapping(region2[0], region2[1], annotation.enums.OverlapType.inside_only)
    assert not annotation_cnv.is_overlapping(region2[0], region2[1], annotation.enums.OverlapType.partial_both)
    assert not annotation_cnv.is_overlapping(region2[0], region2[1], annotation.enums.OverlapType.partial_end)
    assert not annotation_cnv.is_overlapping(region2[0], region2[1], annotation.enums.OverlapType.partial_start)
    assert annotation_cnv.is_overlapping(region2[0], region2[1], annotation.enums.OverlapType.span_whole_only)

    region3 = (12000, 15000)
    assert annotation_cnv.is_overlapping(region3[0], region3[1], annotation.enums.OverlapType.all)
    assert annotation_cnv.is_overlapping(region3[0], region3[1], annotation.enums.OverlapType.inside_only)
    assert not annotation_cnv.is_overlapping(region3[0], region3[1], annotation.enums.OverlapType.partial_both)
    assert not annotation_cnv.is_overlapping(region3[0], region3[1], annotation.enums.OverlapType.partial_end)
    assert not annotation_cnv.is_overlapping(region3[0], region3[1], annotation.enums.OverlapType.partial_start)
    assert not annotation_cnv.is_overlapping(region3[0], region3[1], annotation.enums.OverlapType.span_whole_only)


def test_cnv_get_overlap_with_region(annotation_cnv: CNVRegionAnnotation):
    assert annotation_cnv.get_overlap_with_region(9000, 11000) == 1000
    assert annotation_cnv.get_overlap_with_region(11000, 12000) == 1000
    assert annotation_cnv.get_overlap_with_region(9000, 12000) == 2000
    assert annotation_cnv.get_overlap_with_region(19000, 21000) == 1000
    assert annotation_cnv.get_overlap_with_region(20001, 21000) == 0


def test_cnv_is_position_inside(annotation_cnv: CNVRegionAnnotation):
    assert annotation_cnv.is_position_inside(11000)
    assert not annotation_cnv.is_position_inside(9000)
