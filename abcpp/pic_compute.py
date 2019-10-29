import dendropy
from dendropy.model import continuous
import numpy as np
def pic_compute(tree_sim,Z_singleset,taxa1,order_sim):

    simchar_dict = {}
    keys = ["B.mysticetus", "B.acutorostrata", "C.marginata", "B.borealis",
            "B.physalus", "E.robustus", "B.musculus", "B.omurai",
            "E.australis", "M.novaeangliae", "B.bonaerensis", "B.brydei",
            "B.edeni", "E.glacialis", "E.japonica"]
    for i in range(15):
        simchar_dict[keys[i]] = [Z_singleset[i]]

    simchars = dendropy.ContinuousCharacterMatrix.from_dict(simchar_dict, taxon_namespace=taxa1)
    simpic = continuous.PhylogeneticIndependentConstrasts(tree=tree_sim, char_matrix=simchars)
    sim_ctree = simpic.contrasts_tree(character_index=0,
                                      annotate_pic_statistics=True,
                                      state_values_as_node_labels=False,
                                      corrected_edge_lengths=False)
    sim_pic = []
    sim_label = []
    for nd in sim_ctree.postorder_internal_node_iter():
        sim_pic.append(nd.pic_contrast_standardized)
        sim_label.append(int(nd.label))
    sim_pic_ordered = np.array(sim_pic)[np.argsort(sim_label)]
    return sim_pic_ordered,order_sim