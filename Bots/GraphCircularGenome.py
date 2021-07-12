from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO



def graph_Circular_with_Features_and_Methylation(record, specific_methylation_positions, motif):
    # record = SeqIO.read("../Data/Genome Data/CP035265.gbk", "genbank")        # same result
    gd_diagram = GenomeDiagram.Diagram(record.description)
    gd_track_for_features = gd_diagram.new_track(1, name = "Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()
    gd_methylation_set = gd_track_for_features.new_set()
    for feature in record.features:
        if feature.type != "gene":
            # ignore non-gene features
            continue
        if len(gd_feature_set) % 2 == 0:
            color = colors.blue
        else:
            color = colors.lightblue
        gd_feature_set.add_feature(
            feature = feature,
            color = color,
            label = True,
            label_angle = 45,
        )
        
    for methylation_position in specific_methylation_positions:
        gd_methylation_set.add_feature(
            # feature = format_Motif_to_SeqFeature(methylation_position, motif),
            feature = format_Methlyation_to_SeqFeature(methylation_position),
            color = colors.red,
            lable = True,
            label_angle = 90,
            # label_position = "middle",
            name = "met",
        )

    # Linear example
    gd_diagram.draw(
        format = "linear",
        orientation = "landscape",
        pagesize = "A1",
        fragments = 40,
        start = 0,
        end = len(record),
    )
    gd_diagram.write("MapLinear.pdf","PDF")

    # Circular example
    # gd_diagram.draw(
    #     format = "circular",
    #     circular = True,
    #     pagesize = (200 * cm, 200 * cm), 
    #     start = 0,
    #     end = len(record),
    #     circle_core = 0.9,
    # )
    # gd_diagram.write("MapCircular.pdf", "PDF")


def format_Motif_to_SeqFeature(motif_position, motif):
    return SeqFeature(FeatureLocation(motif_position, motif_position + len(motif)), strand = +1)

def format_Methlyation_to_SeqFeature(methylation_position):
    return SeqFeature(FeatureLocation(methylation_position[0], methylation_position[0] + 1), strand = methylation_position[1])

# record = SeqIO.read("NC_005816.gb", "genbank")
record = SeqIO.read("../Data/Genome Data/CP035265.gbk", "genbank")

gd_diagram = GenomeDiagram.Diagram("TEST PG1 Genome")
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()

for feature in record.features:
    if feature.type != "gene":
        # ignore non-gene features
        continue
    if len(gd_feature_set) % 2 == 0:
        color = colors.blue
    else:
        color = colors.lightblue
    gd_feature_set.add_feature(
        feature = feature,
        color = color,
        label = True,
        label_angle = 45,
    )

# Linear example
gd_diagram.draw(
    format = "linear",
    orientation = "landscape",
    pagesize = "A1",
    fragments = 40,
    start = 0,
    end = len(record),
)
gd_diagram.write("testMapLinear.pdf","PDF")

# Circular example
gd_diagram.draw(
    format = "circular",
    circular = True,
    pagesize = (200 * cm, 200 * cm),
    start = 0,
    end = len(record),
    circle_core = 0.9,
)
gd_diagram.write("testMapCircular.pdf", "PDF")
