import os

for index in range(1, 62):
    if os.path.isfile('/Users/rongpu/Downloads/fringe/plots_compare/original_fringe_{}.png'.format(index)):
        # os.system('mv /Users/rongpu/Downloads/fringe/plots_compare/final_fringe_{}.png /Users/rongpu/Downloads/fringe/plots_compare/{}_final_fringe.png'.format(index, index))
        os.system('mv /Users/rongpu/Downloads/fringe/plots_compare/final_fringe_smooth_{}.png /Users/rongpu/Downloads/fringe/plots_compare/{}_final_fringe_smooth.png'.format(index, index))
        os.system('mv /Users/rongpu/Downloads/fringe/plots_compare/original_fringe_{}.png /Users/rongpu/Downloads/fringe/plots_compare/{}_original_fringe.png'.format(index, index))