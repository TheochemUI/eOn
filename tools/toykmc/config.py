use_sb = True


#tunes barrier size
dEsaddle = .1

#these are the directions that atoms are allowed to move
move_neighbors = [(1,0), (-1,0), (0,1), (0,-1),
        (1,1), (-1,-1), (1,-1), (-1,1)]

#these are the first nearest neighbors for energy purposes
energy_neighbors = [(1,0), (-1,0), (0,1), (0,-1)]

#these are the second nearest neighbors for energy purposes
energy_neighbors_2 = [(1,1), (1,-1), (-1,1), (-1,-1)]
energy_2 = 0.5 #relative strength of 2NN bonds
