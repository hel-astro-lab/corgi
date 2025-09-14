import unittest
import itertools

import pycorgi.twoD as corgi2D
import pycorgitest

class pairwise_moore_communication(unittest.TestCase):

    def test_wout_virtual_tiles(self):
        Nx, Ny = 6, 7
        grid = corgi2D.Grid(Nx, Ny)

        index_space = itertools.product(range(Nx), range(Ny))

        for i, j in index_space:
            t = pycorgitest.MooreTestTile()
            grid.add_tile(t, (i, j))

        grid.pairwise_moore_communication(42)

        def assertTileModes(tile):
            self.assertEqual(len(tile.modes), 8, msg=f"At tile {tile.index}")
            for i, m in enumerate(tile.modes):
                self.assertEqual(m, 42, msg=f"{i}th mode")

        def assertTileCommunications(tile):
            self.assertEqual(len(tile.communications), 8, msg=f"At tile {tile.index}")
            self.assertEqual(len(tile.directions), 8, msg=f"At tile {tile.index}")
            for i, (other_indices, dir) in enumerate(zip(tile.communications, tile.directions)):
                my_indices = tile.index
                expected_other_indices = [(o + d) % N for o, d, N in zip(my_indices, dir, (Nx, Ny))]

                self.assertEqual(other_indices, expected_other_indices, msg=f"{i}th communication")

        def assertPrePostlude(tile):
            self.assertEqual(tile.prelude_mode, 42)
            self.assertEqual(tile.postlude_mode, 42)

        for tile_id in grid.get_local_tiles():
            tile = grid.get_tile(tile_id)
            assertTileModes(tile)
            assertTileCommunications(tile)
            assertPrePostlude(tile)

if __name__ == '__main__':
    unittest.main()

