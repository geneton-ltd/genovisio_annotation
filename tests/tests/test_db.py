import annotation


def test_gnomad():
    annot = annotation.Annotation.load_from_json("tests/data/gnomad.json")

    var_regions = annot.get_common_variability_regions()
    assert len(var_regions) == 1
    assert var_regions[0]["frequency"] == 0.2


def test_genes_annotSV():
    annot = annotation.Annotation.load_from_json("tests/data/genes.json")

    genes = annot.get_annotated_genes()
    assert genes["morbid_genes"] == ["test2"]
    assert genes["associated_with_disease"] == ["test1", "test2"]


def test_get_genes():
    annot = annotation.Annotation.load_from_json("tests/data/genes.json")

    protein_genes = annot.get_genes(gene_type="protein_coding")
    assert len(protein_genes) == 1
    assert protein_genes[0]["gene_name"] == "test3"

    genes = annot.get_genes(overlap=annotation.enums.Overlap.CONTAINED_INSIDE)
    assert len(genes) == 2


def test_count_genes():
    annot = annotation.Annotation.load_from_json("tests/data/genes.json")

    counter = annot.count_gene_types()
    assert counter["protein_coding"] == 1
    assert counter["lncrna"] == 1
    assert counter["pseudogenes"] == 1
    assert counter["other"] == 0
    assert counter["mirna"] == 0
    assert counter["snrna"] == 0
    assert counter["rrna"] == 0


def test_hi_genes():
    annot = annotation.Annotation.load_from_json("tests/data/higenes.json")

    assert len(annot.get_haploinsufficient_genes(annotation.enums.Overlap.ANY, [2, 3])) == 2
    assert len(annot.get_haploinsufficient_genes(annotation.enums.Overlap.ANY, [1, 2, 3])) == 3
    assert len(annot.get_haploinsufficient_genes(annotation.enums.Overlap.ANY, [3])) == 1

    assert len(annot.get_triplosensitivity_genes(annotation.enums.Overlap.ANY, [1, 2, 3])) == 2
    assert len(annot.get_triplosensitivity_genes(annotation.enums.Overlap.ANY, [2, 3])) == 1
    assert len(annot.get_triplosensitivity_genes(annotation.enums.Overlap.ANY, [3])) == 1

    assert annot.get_triplosensitivity_gene_names(annotation.enums.Overlap.ANY, [3]) == ["test2"]


def test_hi_regions():
    annot = annotation.Annotation.load_from_json("tests/data/hiregions.json")

    assert len(annot.get_haploinsufficient_regions(annotation.enums.Overlap.ANY, [2, 3])) == 2
    assert len(annot.get_haploinsufficient_regions(annotation.enums.Overlap.ANY, [1, 2, 3])) == 3
    assert len(annot.get_haploinsufficient_regions(annotation.enums.Overlap.ANY, [3])) == 1

    assert len(annot.get_triplosensitivity_regions(annotation.enums.Overlap.ANY, [1, 2, 3])) == 2
    assert len(annot.get_triplosensitivity_regions(annotation.enums.Overlap.ANY, [2, 3])) == 1
    assert len(annot.get_triplosensitivity_regions(annotation.enums.Overlap.ANY, [3])) == 1

    assert (
        annot.get_triplosensitivity_regions(annotation.enums.Overlap.ANY, [3])[0]["ISCA Region Name"]
        == "test_region_name2"
    )


def test_count_regulatory():
    annot = annotation.Annotation.load_from_json("tests/data/regulatory.json")

    counter = annot.count_regulatory_types()
    assert counter == {
        "enhancer": 1,
        "promoter": 1,
        "open_chromatin_region": 0,
        "CTCF_binding_site": 0,
        "TF_binding_site": 0,
        "curated": 0,
        "flanking_region": 1,
        "silencer": 0,
        "transcriptional_cis_regulatory_region": 0,
        "DNase_I_hypersensitive_site": 0,
        "enhancer_blocking_element": 0,
        "TATA_box": 0,
        "other": 1,
    }
