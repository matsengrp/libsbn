import numpy as np
import sbn

def test_instance():
    inst = sbn.instance('charlie')
    inst.read_newick('data/five_taxon.tre')
    inst.print_status()
    # assert inst.tree_count() == 3
    sbn.f(np.array([3,4]))
    rootsplit_support = inst.rootsplit_support()
    subsplit_support = inst.subsplit_support()

    with open('_build/support.txt', 'w') as fp:
        for support in [rootsplit_support, subsplit_support]:
            support_list = list(support.keys())
            support_list.sort()
            for support in support_list:
                fp.write(support+'\n')
            fp.write('\n')

    inst.read_newick('data/hello.nwk')
    inst.read_fasta('data/hello.fasta')
    inst.beagle_create()
    inst.prepare_beagle_instance()
