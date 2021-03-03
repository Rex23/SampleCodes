from mine import mine

num_row, num_col, num_mine = 16, 16, 128
the_mine = mine(num_row, num_col, num_mine)
the_mine.play()

if __name__ == "__main__":
   print("The mine sweeper game ends.") if not the_mine.success else print("The mine sweeper game ends and you win the game.")