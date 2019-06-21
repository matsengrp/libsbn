import ete3

fp = open('src/doctest.cpp', 'w')

# This tree started life as a tree from the ete documentation:
# t = Tree('((((H,K)D,(F,I)G)B,E)A,((L,(N,Q)O)J,(P,S)M)C);', format=1)

t = ete3.Tree('((((0_1,1_1)1_2,(2_1,3_1)3_2)3_4,4_1)4_5,((5_1,(6_1,7_1)7_2)7_3,(8_1,9_1)9_2)9_5)9_10;', format=1)

preamble = """\
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "driver.hpp"
#include "sbn.hpp"

TEST_CASE("Node") {

Driver driver;

std::vector<std::string> trace;
"""

fp.write(preamble)

fp.write('auto t = driver.ParseString("')
fp.write(t.write(format=9))
fp.write('");\n')

traversal_translator = {
   "preorder": "PreOrder",
   "postorder": "PostOrder"
}

for traversal_type in ["preorder", "postorder"]: # "levelorder"]:
   fp.write("\n// " + traversal_type + ":\n")
   fp.write(f"t->{traversal_translator[traversal_type]}")
   fp.write("([&trace](Node* node) { trace.push_back(node->TagString()); });\n")
   fp.write("CHECK(std::vector<std::string>({")
   fp.write(",".join(['"'+node.name+'"' for node in t.traverse(traversal_type)]))
   fp.write("}) == trace);\n")
   fp.write("trace.clear();\n")

fp.write("\n}\n")

fp.close()
