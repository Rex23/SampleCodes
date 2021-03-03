import random
from collections import deque

class mine:

    def __init__(self, num_row : int, num_col : int, num_mines : int):
        self.num_row = num_row
        self.num_col = num_col
        self.num_mines = num_mines
        self.flags = [[0 for _ in range(num_col)] for _ in range(num_row)] #0: empty, 2: populated, -1: mine, -2: stepped mine
        self.mine_locs = self.generate_mine_locations()
        self.num_total = self.num_row * self.num_col
    
    def play(self):

        counter = self.num_mines

        while counter < self.num_total:
            
            instruction = 0

            input_ = input("click on the location (x, y): ").split(',')

            loc = (int(input_[0]), int(input_[1]))

            if instruction == 0:
                if not self.given_position_instruction(loc, instruction):
                    self.display(1, loc)
                    break
                else:
                    counter += self.mark_empty_cells(loc)
                    self.display(0, loc)
            elif instruction == 1:
                pass

    def display(self, game_over_flag: int, loc: tuple):

        a_map = dict()
        a_map[0] = 'E' if game_over_flag else '*'
        a_map[-1] = 'M' if game_over_flag else '*'
        a_map[-2] = 'M'
        a_map[2] = 'E'
        
        print("The mine after clicking loc: ", loc)
        for row in range(self.num_row):
            print(' '.join([a_map[item] for item in self.flags[row]]))

    def mark_empty_cells(self, loc: tuple) -> int:

        dirs = ((-1, 0), (+1, 0), (0, +1), (0, -1))

        a_q = deque()

        a_q.append(loc)

        self.flags[loc[0]][loc[1]] = 2

        count = 1

        while len(a_q) > 0:

            x, y = a_q.popleft()

            for dir in dirs:
                x_new, y_new = x + dir[0], y + dir[1]

                if x_new >= 0 and x_new < self.num_row and y_new >= 0 and y_new < self.num_col and self.flags[x_new][y_new] == 0:
                    self.flags[x_new][y_new] = 2
                    a_q.append((x_new, y_new))
                    count += 1
        
        return count

    def generate_mine_locations(self) -> set:

        list_mine_locs = random.sample(range(0, self.num_row * self.num_col), self.num_mines)

        locs_2D = set()
        for loc in list_mine_locs:
            x, y = loc // self.num_col, loc % self.num_col
            locs_2D.add( (x, y) )
            self.flags[x][y] = -1
        return locs_2D
    
    def given_position_instruction(self, pos: tuple, instruction: int) -> bool:
        #instruction == 0: sweep it
        #instruction == 1: mark as mine

        if instruction == 0:

            if pos in self.mine_locs:
                self.flags[pos[0]][pos[1]] = -2
                print("Stepped on the mine. Game is over!")
                return False
        
        return True


    