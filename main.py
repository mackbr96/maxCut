import networkx as nx
import cvxopt as cvx
import cvxopt.lapack
import picos as pic
import numpy as np
import matplotlib.pyplot as pylab



def parse(fileName):
    G = nx.Graph()
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
                G.add_edge(int(first), int(second), weight = 1)
    return G


G = parse("ksoso.txt")

N=20

# Generate a G using LCF notation.
#G = nx.LCF_graph(N,[1,3,14],5)


maxcut = pic.Problem()

N = len(G)
#N = 20
X=maxcut.add_variable('X',(N,N),'symmetric')
LL = 1/4.*nx.laplacian_matrix(G).todense()
L=pic.new_param('L',LL)
maxcut.add_constraint(pic.tools.diag_vect(X)==1)

# Constrain X to be positive semidefinite.
maxcut.add_constraint(X>>0)
print("42")
# Set the objective.
maxcut.set_objective('max',L|X)
print("45")
#print(maxcut)

# Solve the problem.
maxcut.solve(verbose = 1,solver='cvxopt')
print("50")
print('bound from the SDP relaxation: {0}'.format(maxcut.obj_value()))
# Use a fixed RNG seed so the result is reproducable.
cvx.setseed(1)
V=X.value
cvxopt.lapack.potrf(V)
for i in range(N):
  for j in range(i+1,N):
    V[i,j]=0

# Do up to 100 projections. Stop if we are within a factor 0.878 of the SDP
# optimal value.
count=0
obj_sdp=maxcut.obj_value()
obj=0
while (count < 1000 or obj < 0.9*obj_sdp):
  r=cvx.normal(N,1)
  x=cvx.matrix(np.sign(V*r))
  o=(x.T*L*x).value[0]
  if o > obj:
    x_cut=x
    obj=o
  count+=1
x=x_cut

# Extract the cut and the seperated node sets.
S1=[n for n in range(N) if x[n]<0]
S2=[n for n in range(N) if x[n]>0]
cut = [(i,j) for (i,j) in G.edges() if x[i]*x[j]<0]

#print( sum(1 for e in cut))
print(sum(G[e[0]][e[1]]['weight'] for e in cut))
print(maxcut.obj_value())


# # Assign colors based on set membership.
# node_colors=[('lightgreen' if n in S1 else 'lightblue') for n in range(N)]

# # Draw the nodes and the edges that are not in the cut.
# nx.draw_networkx(G, pos, node_color=node_colors, edgelist=leave)
# labels={e: '{}'.format(G[e[0]][e[1]]['weight']) for e in leave}
# nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)

# # Draw the edges that are in the cut.
# nx.draw_networkx_edges(G, pos, edgelist=cut, edge_color='r')
# labels={e: '{}'.format(G[e[0]][e[1]]['weight']) for e in cut}
# nx.draw_networkx_edge_labels(G, pos, edge_labels=labels, font_color='r')

# Show the relaxation optimum value and the cut capacity.
rval = maxcut.obj_value()
sval = sum(G[e[0]][e[1]]['weight'] for e in cut)
fig.suptitle(
  'SDP relaxation value: {0:.1f}\nCut value: {1:.1f} = {2:.3f}×{0:.1f}'
  .format(rval, sval, sval/rval), fontsize=16, y=0.97)

# # Show the figure.
# pylab.show()









































print("ok")

