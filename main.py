import networkx as nx
import cvxopt as cvx
import cvxopt.lapack
import picos as pic


def parse(fileName):
    graph = nx.Graph()
    with open(fileName) as f:
        lines = f.readlines()
        for line in lines:
            if "graph{" not in line and "}" not in line:
                line = line.replace(" ", "")
                line = line.replace(";", "")
                line = line.strip()
                first = ""
                second = ""
                i = 0
                while line[i] != "-":
                    first += line[i]
                    i += 1
                i += 2
                while i < len(line):
                    second += line[i]
                    i += 1
                graph.add_edge(int(first), int(second))
    return graph

                    
graph = parse("ksoso.txt")

maxcut = pic.Problem()

N = len(graph)
X=maxcut.add_variable('X',(N,N),'symmetric')
LL = 1/4.*nx.laplacian_matrix(graph).todense()









































print("ok")

