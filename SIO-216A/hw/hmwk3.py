""" Methods and the like related to Homework 3"""
import numpy as np

from IPython import embed

class Ants():
    def __init__(self, nants, N, M, nest,
                 trace_food=1e-3, trace_nest=1e-3,
                 sniff_food=0.3, sniff_nest=0.3):
        # Grid
        self.N = N
        self.M = M
        # Pheremones
        self.trace_food = trace_food
        self.trace_nest = trace_nest
        self.sniff_food = sniff_food
        self.sniff_nest = sniff_nest
        #
        self.nants = nants
        self.x = nest[0] * np.ones(nants, dtype=int)
        self.y = nest[1] * np.ones(nants, dtype=int)
        self.goal = np.ones(nants, dtype=bool)
        #
        self.y_neighbor = [1, 1, 1, 0, 0, -1, -1, -1]
        self.x_neighbor = [-1, 0, 1, -1, 1, -1, 0, 1]
        self.theta_neighbor = np.arctan2(-1*np.array(self.y_neighbor),
                                         -1*np.array(self.x_neighbor))

    def update_xy(self, theta, subset=None):
        if subset is None:
            subset = np.arange(len(theta))
        # Ugly loop on ants
        for ii, iant in enumerate(subset):
            idx_th = int(np.argmin(np.abs(theta[ii] - self.theta_neighbor)))
            #
            newx = self.x[iant] + self.x_neighbor[idx_th]
            newy = self.y[iant] + self.y_neighbor[idx_th]
            #
            if newx == self.N:
                newx = self.N - 2
            elif newx == 0:
                newx = 1
            if newy == self.M:
                newy = self.M - 2
            elif newy == 0:
                newy = 1
            # Assign
            self.x[iant] = newx
            self.y[iant] = newy

    def grid_ants(self, atype):
        grid = np.zeros((self.N, self.M), dtype=int)
        if atype == 'food':
            for x,y in zip(self.x[self.goal], self.y[self.goal]):
                grid[x, y] += 1
        elif atype == 'nest':
            for x,y in zip(self.x[~self.goal], self.y[~self.goal]):
                grid[x, y] += 1
        #
        return grid


def decay_phenom(p_food, p_nest, p_decay=2.):
    p_food /= p_decay
    p_nest /= p_decay

def step_to_food(rstate, ants, food, p_food,
                        random_step=0.1):
    """

    Parameters
    ----------
    rstate
    ants : Ants
    food
    p_food
    random_step

    Returns
    -------

    """
    # Cut on to food
    food_seeking = ants.goal
    # Go straight towards food?
    to_food = (rstate.rand(ants.nants) < ants.sniff_food) & food_seeking
    if np.any(to_food):
        theta = np.arctan2((ants.y[to_food]-food[1]),(ants.x[to_food]-food[0]))
        # Choose new cell
        ants.update_xy(theta, subset=np.where(to_food)[0])

    # Go random?
    go_random = (rstate.rand(ants.nants) < random_step) & np.invert(to_food) & food_seeking
    if np.any(go_random):
        theta = 2*np.pi*rstate.rand(np.sum(go_random)) - np.pi
        # Choose new cell
        ants.update_xy(theta, subset=np.where(go_random)[0])

    # Use pheremones
    use_pherome = np.invert(go_random) & np.invert(to_food) & food_seeking
    for ii in np.where(use_pherome)[0]:
        vals = []
        for ix, iy in zip(ants.x_neighbor, ants.y_neighbor):
            vals.append(p_food[ants.x[ii]+ix, ants.y[ii]+iy])
        ival = np.argmax(vals)
        ants.x[ii] += ants.x_neighbor[ival]
        ants.y[ii] += ants.y_neighbor[ival]

    # Add pheremone
    for x, y in zip(ants.x, ants.y):
        p_food[x,y] += ants.trace_food

    # Update those that got there!
    arrived = (ants.x == food[0]) & (ants.y == food[1])
    if np.any(arrived):
        ants.goal[arrived] = False


def step_to_nest(rstate, ants, nest, p_nest, random_step=0.1):
    """

    Parameters
    ----------
    rstate
    ants : Ants
    nest
    p_nest
    random_step

    Returns
    -------

    """
    # Cut on to nest
    nest_seeking = ~ants.goal
    # Go straight towards nest?
    to_nest = (rstate.rand(ants.nants) < ants.sniff_nest) & nest_seeking
    if np.any(to_nest):
        theta = np.arctan2((ants.y[to_nest] - nest[1]), (ants.x[to_nest] - nest[0]))
        # Choose new cell
        ants.update_xy(theta, subset=np.where(to_nest)[0])

    # Go random?
    go_random = (rstate.rand(ants.nants) < random_step) & np.invert(to_nest) & nest_seeking
    if np.any(go_random):
        theta = 2 * np.pi * rstate.rand(np.sum(go_random)) - np.pi
        # Choose new cell
        ants.update_xy(theta, subset=np.where(go_random)[0])

    # Use pheremones
    use_pherome = np.invert(go_random) & np.invert(to_nest) & nest_seeking
    for ii in np.where(use_pherome)[0]:
        vals = []
        for ix, iy in zip(ants.x_neighbor, ants.y_neighbor):
            vals.append(p_nest[ants.x[ii] + ix, ants.y[ii] + iy])
        ival = np.argmax(vals)
        ants.x[ii] += ants.x_neighbor[ival]
        ants.y[ii] += ants.y_neighbor[ival]

    # Add pheremone
    for x, y in zip(ants.x, ants.y):
        p_nest[x, y] += ants.trace_nest

    # Update those that got there!
    arrived = (ants.x == nest[0]) & (ants.y == nest[1])
    if np.any(arrived):
        ants.goal[arrived] = True


def run_forward(rstate, nsteps, ants, food, nest, p_food, p_nest,
                random_step=0.1):
    for ii in range(nsteps):
        # Take steps to food
        step_to_food(rstate, ants, food, p_food, random_step=random_step)
        # Take steps to nest
        step_to_nest(rstate, ants, nest, p_nest, random_step=random_step)
        # Decay
        decay_phenom(p_food, p_nest)

if __name__ == '__main__':
    N, M = 100, 100
    nest = (10, 10)
    food = (35, 75)
    p_food = np.zeros((N, M))
    rstate = np.random.RandomState()
    #
    ants = Ants(100, N, M, nest)
    for ii in range(500):
        random_step_to_food(rstate, ants, food, p_food,
                                  food_sniff=0.7, random_step=0.1)
