""" Magic Golden Coin problem """
import numpy as np

from IPython import embed

def rowcols_grid(grid, axis):
    trufalse = ((np.sum(grid, axis=axis)) % 2) == 1
    return trufalse.astype(int)

def nrowcol_to_rowcolpos(nrowcol, ngrid=2):
    npos = int(np.log2(ngrid))

    binrow = bin(nrowcol[1])[2:]
    binrow = '0'*(npos-len(binrow)) + binrow
    bincol = bin(nrowcol[0])[2:]
    bincol = '0'*(npos-len(bincol)) + bincol

    colpos = np.array([int(item) for item in bincol])
    rowpos = np.array([int(item) for item in binrow])
    # These read left to right
    return colpos, rowpos

def rowcolpos_to_nrowcol(colpos, rowpos):
    twos = [2**kk for kk in range(len(colpos))]
    twos.reverse()
    nrow = np.sum([two*item for two,item in zip(twos,rowpos.tolist())])
    ncol = np.sum([two*item for two,item in zip(twos,colpos.tolist())])
    return ncol, nrow

def flip(cell):
    if cell == 1:
        return 0
    else:
        return 1

def encode(grid, magic_pos, ngrid=2, debug=False):
    magic_colpos, magic_rowpos = nrowcol_to_rowcolpos(
        magic_pos,  ngrid=ngrid)
    ogrid = grid.copy()

    # Check if we are correct already
    ogrid_col, ogrid_row = decode(ogrid, ngrid=ngrid)
    if np.all(magic_colpos == ogrid_col) & np.all(magic_rowpos == ogrid_row):
        return grid

    # Brute force me
    flip_row = None
    for row in range(ngrid):
        grid = ogrid.copy()
        grid[row,0] = flip(ogrid[row,0])
        grid_col, grid_row = decode(grid, ngrid=ngrid)
        print('row:', grid_row)
        if np.all(magic_rowpos == grid_row):
            flip_row = row
            break

    flip_col = None
    for col in range(ngrid):
        grid = ogrid.copy()
        grid[0, col] = flip(ogrid[0, col])
        grid_col, grid_row = decode(grid, ngrid=ngrid)
        print('col:', grid_col)
        if np.all(magic_colpos == grid_col):
            flip_col = col
            break

    # Do it
    grid = ogrid.copy()
    try:
        grid[flip_row,flip_col] = flip(ogrid[flip_row,flip_col])
    except ValueError:
        embed(header='70 of coin')

    if debug:
        embed(header='debug: 73 of coin')

    # Return
    return grid

def decode(grid, ngrid=2):
    npos = int(np.log2(ngrid))
    transforms = np.arange(ngrid+1).tolist()
    transforms.reverse()

    drop_list = [0] + [2**item for item in range(npos+1)]
    for drop in drop_list:
        transforms.remove(drop)
    #embed(header='87 of coin')
    #
    all_col_pos = rowcols_grid(grid, 0)
    all_row_pos = rowcols_grid(grid, 1)
    col_pos = all_col_pos[-npos:]
    row_pos = all_row_pos[0:npos]

    # Transform
    for tt in range(ngrid-npos-1):
        bint = bin(transforms[tt])[2:]
        bint = '0'*(npos-len(bint)) + bint
        #if tt > 0:
        #    embed(header='91 of transform')
        #
        col_trans = all_col_pos[-npos-1-tt]
        if col_trans == 1:
            for kk, tval in enumerate(bint):
                if tval == '1':
                    col_pos[kk] = flip(col_pos[kk])
        # Rows
        row_trans = all_row_pos[npos+tt]
        if row_trans == 1:
            for kk, tval in enumerate(bint):
                if tval == '1':
                    row_pos[kk] = flip(row_pos[kk])
    '''
    if ngrid > 2:
        # Columns
        col_flip_all = all_col_pos[-npos-1]
        if col_flip_all == 1:
            for kk in range(col_pos.size):
                col_pos[kk] = flip(col_pos[kk])
        # Rows
        row_flip_all = all_row_pos[npos]
        if row_flip_all == 1:
            for kk in range(row_pos.size):
                row_pos[kk] = flip(row_pos[kk])
    '''
        
    # Return
    return col_pos, row_pos


if __name__ == '__main__':
    # 4x4
    '''
    ngrid = 4
    rand_grid = np.random.randint(0, high=2, size=(ngrid, ngrid))
    col, row = decode(rand_grid, ngrid=ngrid)
    '''

    # Random trials
    nrand = 50
    ngrid = 64
    for ss in range(nrand):
        # Grid
        rand_grid = np.random.randint(0, high=2, size=(ngrid, ngrid))
        # Coin
        rand_coin = np.random.randint(0, high=ngrid, size=2).tolist()
        
        # Encode
        new_grid = encode(rand_grid, rand_coin, ngrid=ngrid,
                          debug=False)

        # Decode
        col, row = decode(new_grid, ngrid=ngrid)
        ncol, nrow = rowcolpos_to_nrowcol(col, row)

        #
        try:
            assert [ncol, nrow] == rand_coin
        except:
            embed(header='124 of coin')
        print(f"passed: ss={ss}")

    '''
    grid = np.array( [[1,0], [0,1]] ) # col=0, row=0
    magic_pos = [0,0]
    new_grid = encode(grid, magic_pos)

    col, row = decode(new_grid)
    ncol, nrow = rowcolpos_to_nrowcol(col, row)

    print(f'Magic coin at col={ncol}, row={nrow}; Ans={magic_pos}')
    embed(header='75 of coin')

    # Test magic pos
    magic_rowcolpos = nrowcol_to_rowcolpos([2,3], ngrid=4)


    # Test translator
    grid = np.array( [[0,1], [0,1]] ) # col=0, row=1

    grid = np.array( [[1,1], [0,1]] ) # col=0, row=0

    #grid = np.array( [[0,0], [1,1]] )

    #grid = np.array( [[0,1], [1,1]] )  # col=0, row=1
    #grid = np.array( [[0,1], [0,1]] )


    col, row = decode(grid)
    print(f'col={col}, row={row}')

    embed(header='25 of coin')
    '''