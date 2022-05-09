set term pdf

unset g
set view map
set pm3d at b map
set dgrid3d 200,100,2

set o "results/boundary.pdf"
file = "results/boundary.txt"
splot file u 1:2:3 not

set o "results/phi.pdf"
file = "results/phi.txt"
splot file u 1:2:3 not

set o "results/dissociation.pdf"
file = "results/dissociation.txt"
splot file u 1:2:3 not
