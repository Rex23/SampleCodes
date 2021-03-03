from mine import mine

num_row, num_col, num_mine = 5, 3, 5
the_mine = mine(num_row, num_col, num_mine)
the_mine.play()

if __name__ == "__main__":
   print("The mine sweeper game ends.")