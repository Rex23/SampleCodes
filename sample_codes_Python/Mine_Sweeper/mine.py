import random
import copy
from collections import deque

class mine:

    def __init__(self, num_row : int, num_col : int, num_mines : int):
        self.num_row = num_row
        self.num_col = num_col
        self.num_mines = num_mines
        self.flags = [[0 for _ in range(num_col)] for _ in range(num_row)] # >=0: empty, -1: mine, -2: stepped mine
        self.populated = [[0 for _ in range(num_col)] for _ in range(num_row)] # 0: no, 1: yes
        self.mine_locs = self.generate_mine_locations()
        self.generate_flags()
        self.flags_ori = copy.deepcopy(self.flags)
        self.num_total = self.num_row * self.num_col
        self.success = True
    
    def generate_flags(self):
        
        dirs = ((-1, 0), (+1, 0), (0, +1), (0, -1), (-1, -1), (-1, +1), (+1, -1), (+1, +1))
        
        for i in range(self.num_row):
            for j in range(self.num_col):
                if self.flags[i][j] != -1:
                    count = 0
                    for dir in dirs:
                        i_new, j_new = i + dir[0], j + dir[1]
                        if i_new >= 0 and i_new < self.num_row and j_new >= 0 and j_new < self.num_col and self.flags[i_new][j_new] == -1:
                            count += 1
                    self.flags[i][j] = count
    
    def play(self):

        counter = self.num_mines

        while counter < self.num_total:
            
            input_ = input("click on the location x, y, instruction: ").split(',')

            loc = (int(input_[0]), int(input_[1]))

            instruction = input_[2].strip()

            if instruction == "R":
                if not self.given_position_instruction(loc, instruction):
                    self.success = False
                    self.display(1, loc)
                    break
                else:
                    counter += self.mark_empty_cells(loc)
                    self.display(0, loc)
            elif instruction == "F":
                self.mark_empty_cell_with_red_flag(loc)
                self.display(0, loc)
    
    def mark_empty_cell_with_red_flag(self, loc: tuple):

        if self.flags[loc[0]][loc[1]] != -2:
            self.flags[loc[0]][loc[1]] = -2
            self.populated[loc[0]][loc[1]] = 1
        else:
            self.flags[loc[0]][loc[1]] = self.flags_ori[loc[0]][loc[1]]
            self.populated[loc[0]][loc[1]] = 0

    def display(self, game_over_flag: int, loc: tuple):
        
        print("The mine after clicking loc: ", loc)
        for i in range(self.num_row):
            for j in range(self.num_col):
                
                if self.flags[i][j] == -1:
                    if game_over_flag:
                        print("M", end = ", ") if j != self.num_col - 1 else print("M")
                    else:
                        print("*", end = ", ") if j != self.num_col - 1 else print("*")
                elif self.flags[i][j] != -2:
                    if game_over_flag:
                        print(self.flags[i][j], end = ", ") if j != self.num_col - 1 else print(self.flags[i][j])
                    else:
                        if self.populated[i][j] == 1:
                            print(self.flags[i][j], end = ", ") if j != self.num_col - 1 else print(self.flags[i][j])
                        else:
                            print("*", end = ", ") if j != self.num_col - 1 else print("*")
                elif self.flags[i][j] == -2:
                    if game_over_flag:
                        if self.flags_ori[i][j] == -1:
                            print("M", end = ", ") if j != self.num_col - 1 else print("M")
                        else:
                            print(self.flags_ori[i][j], end = ", ") if j != self.num_col - 1 else print(self.flags_ori[i][j])
                    else:
                        print("F", end = ", ") if j != self.num_col - 1 else print("F")

    def mark_empty_cells(self, loc: tuple) -> int:

        dirs = ((-1, 0), (+1, 0), (0, +1), (0, -1), (-1, -1), (-1, +1), (+1, -1), (+1, +1))

        a_q = deque()

        a_q.append(loc)

        self.populated[loc[0]][loc[1]] = 1

        count = 1

        while len(a_q) > 0:

            x, y = a_q.popleft()

            for dir in dirs:
                x_new, y_new = x + dir[0], y + dir[1]

                if x_new >= 0 and x_new < self.num_row and y_new >= 0 and y_new < self.num_col \
                and self.flags[x_new][y_new] != -1 \
                and self.populated[x_new][y_new] == 0:
                    self.populated[x_new][y_new] = 1
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
    
    def given_position_instruction(self, pos: tuple, instruction: str) -> bool:

        if instruction == "R":

            if pos in self.mine_locs:
                print("Stepped on the mine. Game is over!")
                return False
        
        return True


    