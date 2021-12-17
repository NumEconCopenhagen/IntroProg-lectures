import numpy as np 
import random  
import math
import types
import seaborn 
import matplotlib.pylab as plt
import matplotlib.pyplot as pyplot
import matplotlib.image as img
import imageio
import os
import shutil



class Simulation:

    def __init__(self):

        # a. Baseline parameters
        self.n_dim = 10 
        self.check_prints = False
        self.cutoff = 0.7
        self.T = 5
        self.t = 0
        self.share_pop = 0.8
        self.share_A = 0.5
        self._set_dimensions() # Compute dimensions of landscape and population 
        
        # b. Containers for agents and squares
        self.citizens = []
        self.squares = []
        self.figs = []
        self.is_setup = False
    
    def get_model_params(self):
        ''' Obtain the current parameter settings of the simulation

        Returns:
            mp (SimpleNamespace): an object with the same attributes as the simulation parameters
        ''' 

        mp = types.SimpleNamespace()
        
        # a. Convert the existing parameters to a dict
        all_attributes = self.__dict__ 
        
        # b. Loop over attributes obtained from simulation (self) and assign them to mp 
        for param,value in all_attributes.items():
            if isinstance(value, (list, set)):
                continue # Exclude the lists of all citizens and squares
            if param not in ['citiznes', 'squares', 'figs', 'is_setup', '_set_dimensions', 't', 'n_A', 'n_B', 'n_squares', 'n_pop']: 
                setattr(mp, param, value)
        return mp

    def setup_sim(self, mp):
        ''' Setup the simulation by updating settings, create landscape and citizens

        The lists of squares (the landscape) and citizens created here are stored as attributes in the simulation object.

        Args:
            mp (SimpleNamespace): an object containing all model parameters. Obtained by get_model_params(). 
        '''

        
        # a. Update settings according to input mp
        if self.is_setup:
            self.citizens = []
            self.squares = []
            self.figs = []
            self.t = 0
            
        mp_dict = mp.__dict__
        for param,value in mp_dict.items():
            setattr(self, param, value)
        self._set_dimensions()

        # b. Create the squares of the landscape
        for y in range(self.n_dim):
            for x in range(self.n_dim):
                sq = Square(x_coord = x, y_coord = y)
                self.squares.append(sq)

        # b.i Add to each square a list of it's neighbouring squares as an attribute 
        for sq in self.squares:
            sq.neighbours = self.get_neighbours(sq)

        # c. Create the citizens of the simulation
        free_squares = list(range(self.n_squares))  # All squares are free at the start of the simulation
        random.shuffle(free_squares)

        for i in range(self.n_pop):
            # c.i Create a citizen and assign to it an initial square
            citizen = Citizen()
            citizen.square = self.squares[free_squares.pop()]

            # c.ii Define the type of citizen i
            if i < self.n_A :
                citizen.type = 'A'
            else :
                citizen.type = 'B'

            # c.iii Update status of the square to citizen i
            citizen.square.is_free = False 
            citizen.square.citizen = citizen
            self._modify_neighbours_status(citizen, citizen.square, 1)

            # c.iv Add citizen to list of all citizens 
            self.citizens.append(citizen)
        
        # d. Create figure path and record initial state of simulation 
        cwd = os.getcwd()
        self.fig_dir = os.path.join(cwd, 'figs')
        if os.path.exists(self.fig_dir):
            shutil.rmtree(self.fig_dir)
        os.mkdir(self.fig_dir)
        self.plot_state(save_plot=True, suppress_plot=True)
        
        self.is_setup = True

    def run_sim(self):
        ''' Run the simulation over T periods. 
        '''

        # The main loop over the simulation 
        for t in range(1,self.T+1):
            self.t = t

            # a. Get all available squares given citizen type. Note we are creating sets by set comprehension.
            free_squares_A = {sq for sq in self.squares if sq.is_free 
                              & self.is_citizen_content('A', sq)}
            free_squares_B = {sq for sq in self.squares if sq.is_free
                              & self.is_citizen_content('B', sq)}

            free_squares = {'A': free_squares_A, 'B': free_squares_B}
            random.shuffle(self.citizens) # Randomize list of citizens every period to avoid systematic 

            # b. Get all discontented citizens at the beginning of iteration t
            discontented = []
            for citizen in self.citizens:
                content = self.is_citizen_content(citizen.type, citizen.square)   
                if not content:
                    discontented.append(citizen)

            # c. Let all discontented citizens move
            for dc in discontented:                  
                # c.i Pop a free square (if any exists) from the relevant set 
                # and remove it from the other type's available squares
                if free_squares[dc.type]:
                    new_square = free_squares[dc.type].pop()
                    other_type = self._get_other_type(dc.type)
                    free_squares[other_type].discard(new_square)

                    # c.ii Move citizen dc away from its current square by deleting dc's presence there.  
                    self._modify_neighbours_status(dc, dc.square, -1)
                    dc.square.citizen = None
                    dc.square.is_free = True 

                    # c.iii Move dc to new square by setting dc's square-attribute to new_square.  
                    dc.square = new_square 
                    dc.square.is_free = False 
                    dc.square.citizen = dc
                    self._modify_neighbours_status(dc, dc.square, 1)

            # d. Save a plot of the simulation state for later GIF-creation
            self.plot_state(save_plot=True, suppress_plot=True)

            if self.check_prints:
                print(f'Iteration {t}: discontented citizens = {len(discontented)}')

    def is_citizen_content(self, citizen_type, square):
        ''' Test if a citizen is content with a given square based on the square's neighbours.
        A citizen is content if the share of neighbors with similar type is above the cut-off level.
        Note: a citizen will always be content with a square that has 0 neighbouring citizens.

        Args: 
            citizen_type (str): the type of citizen in question
            square (Square): the square that the citizen considers
        
        Returns:
            (bool): True if the citizen is content with the square.
        '''

        # a. Create the attribute label for the count of similar type neighbors at the square
        same_types = getattr(square, citizen_type + '_neighbours')
        n_neighbours = square.A_neighbours + square.B_neighbours

        # b. Test if neighbour-share is above cut-off 
        if n_neighbours == 0:
            return True 
        else:
            return (same_types/n_neighbours) >= self.cutoff

    def get_idx_1d(self, x, y):
        '''
        Convert x-y coordinates to a 1D index
        '''
        return y*self.n_dim + x

    def get_neighbours(self, square):
        ''' Create a list of neighboring squares to a given square

        Args:
            square (Square): a square in the landscape of the model
        
        Returns:
            neighbours (list): all the adjacent squares to the input square
        '''

        # a. Create lists of possible neighbours to a given square in each direction
        # a.i Set of adjancent coordinates along the y-axis
        nr = 1
        xc = square.x_coord
        yc = square.y_coord
        
        diam = 2*nr + 1
        x_coords = []
        y_coords = []
        for r in range(diam):
            xvec = list(range(xc-nr, xc+nr+1))
            yvec = [yc - (nr-r)]*diam
            # Exclude the middle square - not a neighbour
            if r == nr:
                del xvec[nr]
                del yvec[nr]
            x_coords = x_coords + xvec
            y_coords = y_coords + yvec
    
        neighbours = []

        # b. Loop through possible neighbours. Discard coordinates that are outside the board and append the rest.
        for x, y in zip(x_coords, y_coords):
            if (x < 0) | (y < 0) | (x >= self.n_dim) | (y >= self.n_dim):
                continue
            else:
                idx = self.get_idx_1d(x, y)
                neighbours.append(self.squares[idx])

        return neighbours

    def plot_state(self, show_t = None, save_plot=False, suppress_plot=False):
        ''' Plot the current state of the simulation

        Args:
            show_t (int): display state of simulation at iteration t
            save_plot (bool): if true, the created plot is stored in figs/
            suppress_plot (bool): if true, the plot is not displayed
        '''

        # a. Display a saved figure if show_t is provided
        if show_t is not None:
            if show_t >= 0 & show_t < len(self.figs):
                ipth = self.figs[show_t]
                image = img.imread(ipth)
                _, ax = pyplot.subplots(figsize=(9,11))
                ax.imshow(image)
                ax.axis('off') 
                return

        # b. Create a 2d numpy array with square states 
        data = np.empty((self.n_dim, self.n_dim))

        for s in self.squares:
            if s.is_free:
                data[s.x_coord, s.y_coord] = 0
            elif s.citizen.type == 'A':
                data[s.x_coord, s.y_coord] = 1
            elif s.citizen.type == 'B':
                data[s.x_coord, s.y_coord] = 2
            else:
                raise Exception('plot_state: Unknown agent type at square')

        # c. Use a Seaborn heatmap for plotting the simulation landscape
        ax = seaborn.heatmap(data, linewidth=0, xticklabels=False,
                             yticklabels=False, cbar=False, cmap='YlGnBu')
        ax.set_title(f'Landscape at iteration {self.t}')
        
        # d. Save plot for creating a GIF post simulation
        if save_plot:
            fig_name = f'fig_state_{self.t}.png'
            fig_name = os.path.join(self.fig_dir, fig_name)
            self.figs.append(fig_name)
            plt.savefig(fig_name)
        
        if suppress_plot:
            plt.close()
        else:
            plt.show()
    
    def make_gif(self, clean_up=False):
        ''' Create a GIF of the recorded simulation landscapes
        '''

        # a. Collect recorded images into a GIF
        with imageio.get_writer('tmp.gif', mode='I') as writer:
            for fig in self.figs:
                image = imageio.imread(fig)
                writer.append_data(image)

        # b. Adjust frames pr second
        gif = imageio.mimread('tmp.gif')
        imageio.mimsave('simulation.gif', gif, fps=3)

        # c. Clean up
        if clean_up:
            self.figs.append('tmp.gif')
            for fn in self.figs:
                os.remove(fn)
            self.figs = []


    def print_neighbours(self, square):
        ''' Print out the state of neighbouring squares to a given square

        Args:
            square (Square): the square for which a print out will be made
        '''

        sq_type = square.citizen.type if square.citizen is not None else 'None'
        print(f'Citizen at ({square.x_coord},{square.y_coord}) is {sq_type}')
        for neighbour in square.neighbours:
            ntype = neighbour.citizen.type if neighbour.citizen is not None else 'None'
            print(
                f'Neighbor at ({neighbour.x_coord},{neighbour.y_coord}) is {ntype}')

    # ========================== 
    # Internal helper methods 
    # ==========================

    def _modify_neighbours_status(self, citizen, square, delta):
        ''' Change the number of neighbours of a given type that a square has.
        Loops through the neighbours of the input square. Adds one (subtracts one) to their count
        of neighbours that has type 'citizen.type' when citizen moves to (moves from) the square.

        Args:
            citizen (Citizen): the agent that is moving in or away from a square.
            square (Square): the square that citizen is moving to or from.     
            delta (int): +1 if citizen is moving to a square. -1 if moving away from square.
        '''

        # Define the type of neighbours to be updated, A_neighbours or B_neighbours, given by citizen type
        my_type = citizen.type + '_neighbours'

        # + 1 or - 1 (delta = -1) to neighbors' counts depending on moving away or moving in.
        for neighbour in square.neighbours:
            n_of_my_type = delta + getattr(neighbour, my_type)
            setattr(neighbour, my_type, n_of_my_type)

    def _get_other_type(self, citizen_type):
        ''' Get the opposite type of input citizen_type
        '''
        other_type = 'A' if citizen_type == 'B' else 'B'
        return other_type

    def _set_dimensions(self):
        ''' Calculate the dimensions of the landscape and population. Internal function.
        '''
        self.share_B = 1 - self.share_A
        self.n_squares = self.n_dim * self.n_dim
        self.n_pop = round(self.n_squares * self.share_pop)
        self.n_A = round(self.n_pop * self.share_A)
        self.n_B = self.n_pop - self.n_A

class Square:
    ''' Objects that make up the landscape of the simulation. 
    Is defined by an (x,y) coordinate. Holds references to neighbouring squares. If inhabitated, squares hold a reference to the citizen living there. 
    '''

    def __init__(self, x_coord = np.nan, y_coord = np.nan):
        ''' A square is initiated with a (x,y) coordinate
        '''

        # a. Initiate properties of a square
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.neighbours = [] 
        self.is_free = True 
        self.citizen = None
        self.A_neighbours = 0
        self.B_neighbours = 0
    

class Citizen:
    ''' The agents of the simulation. 
    Note: all behavior is put into the Simulation class
    '''

    def __init__(self):

        # a. Attributes of citizens
        self.type = None 
        self.square = None



