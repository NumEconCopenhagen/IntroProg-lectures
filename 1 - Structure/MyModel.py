import solve
import simulate

class MyModelClass:

    def __init__(self,x=1):
        
        self.x = x
        self.y = None

        print('setup done')

    solve = solve.all
    simulate = simulate.all

    def save_results(self):

        with open('results/xy.txt', 'w') as file_ref: # 'w' is for 'write'
            file_ref.write(f'x = {self.x}\n')
            file_ref.write(f'y = {self.y}\n')