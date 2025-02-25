from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

##########################################################################################

skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8.fits')

expnum_list = np.array([240789, 241069, 241188, 242471, 243261, 243665, 247935, 251549, 252286, 252332, 277324, 279452, 338730, 349651, 349769, 351226, 351743, 354990, 355400, 360689, 361509, 362059, 362437, 362450, 362944, 364693, 367484, 369471, 370222, 370674, 370674, 376947, 379276, 379901, 380009, 380009, 381919, 386690, 386692, 386692, 387710, 387740, 388113, 388177, 388177, 388596, 388596, 389485, 389485, 390171, 390171, 390252, 390500, 390500, 390594, 390894, 390894, 390989, 391381, 391728, 392465, 392754, 392755, 392778, 392816, 397487, 400429, 400445, 400799, 401940, 402276, 407636, 407965, 409455, 425811, 425857, 426232, 426287, 426287, 432000, 432827, 433179, 433308, 449012, 449607, 449625, 449633, 463545, 463546, 466040, 473885, 481386, 482687, 483663, 483738, 484576, 485415, 486818, 486851, 488787, 493258, 494478, 495049, 502695, 502706, 509522, 511458, 511498, 511521, 514512, 514844, 521619, 521621, 522043, 534040, 535123, 535130, 535227, 547165, 547254, 548280, 553662, 563283, 572207, 572690, 573547, 576442, 576539, 576583, 577017, 577711, 578271, 578595, 579151, 579468, 579519, 579768, 579937, 581277, 585951, 586954, 587126, 587129, 587350, 588650, 588950, 589901, 590219, 591555, 595805, 595867, 595875, 595884, 597313, 598659, 599759, 600174, 600965, 603027, 603737, 604733, 605396, 605624, 611909, 618342, 621204, 624616, 624959, 625350, 625787, 626051, 626487, 630895, 634219, 635653, 635724, 636002, 636004, 646583, 646591, 649500, 650128, 659510, 660499, 660512, 660520, 660574, 660581, 661579, 662403, 663769, 674055, 674764, 677312, 677742, 678353, 679871, 680146, 680166, 680182, 682051, 682111, 686483, 690551, 692803, 695627, 695658, 696578, 698937, 699225, 703841, 705147, 705155, 712799, 720061, 720505, 720825, 721180, 721181, 721198, 721214, 721221, 721232, 721534, 721555, 723111, 723813, 723815, 724166, 724177, 724586, 724894, 725262, 725286, 725288, 730118, 730395, 731076, 731218, 731263, 731269, 731469, 731484, 731488, 731490, 731543, 731919, 745744, 747069, 762941, 763211, 763714, 764237, 769549, 769676, 769984, 770589, 773418, 780162, 780762, 780962, 781378, 781434, 782831, 783610, 783668, 783933, 783967, 784302, 789414, 789906, 790230, 790982, 792757, 793040, 793554, 793581, 793835, 801602, 801966, 802999, 808616, 819582, 819856, 820175, 820199, 829678, 830068])

np.random.seed(1674)
expnum_list = np.random.choice(expnum_list, size=50, replace=False)
print(expnum_list)

f = open("/global/cfs/cdirs/desi/users/rongpu/plots/dr9dev/pupil_pattern/new_images/compare_new_images.html", "w")
f.write('<html>\n')
f.write('<table>\n')
f.write('<th>New</th>\n')
f.write('<th>New + sky_corerction</th> \n')
f.write('<th>Original</th>\n')
f.write('<th>Original + sky_corerction</th>\n')

for expnum in expnum_list:

    skyrun_index = np.where(skyrun['expnum']==expnum)[0][0]
    band = skyrun['filter'][skyrun_index]
    run = skyrun['run'][skyrun_index]
    plot_fn1 = '{}_{}_image_{}_new.png'.format(band, run, expnum)
    plot_fn2 = '{}_{}_image_{}_new_sky_corrected.png'.format(band, run, expnum)
    plot_fn3 = '{}_{}_image_{}_original.png'.format(band, run, expnum)
    plot_fn4 = '{}_{}_image_{}_original_sky_corrected.png'.format(band, run, expnum)
    
    f.write('<tr>\n')

    f.write('<td><a href=\'{}\'><img src=\'{}\' width=\'400\'></a></td>\n'.format(plot_fn1, plot_fn1))
    f.write('<td><a href=\'{}\'><img src=\'{}\' width=\'400\'></a></td>\n'.format(plot_fn2, plot_fn2))
    f.write('<td><a href=\'{}\'><img src=\'{}\' width=\'400\'></a></td>\n'.format(plot_fn3, plot_fn3))
    f.write('<td><a href=\'{}\'><img src=\'{}\' width=\'400\'></a></td>\n'.format(plot_fn4, plot_fn4))
    
    # f.write('<td><a href=\'lg_regular_cp-sky_corrected/c4d_180221_072207_ooi_g_ls9.png\'><img src=\'lg_regular_cp-sky_corrected/c4d_180221_072207_ooi_g_ls9.png\' width=\'400\'></a></td>\n')
    
    f.write('</tr>\n')
        
f.write('</table>\n')
f.close()
