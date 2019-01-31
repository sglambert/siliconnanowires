import os
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import tb
import matplotlib.pyplot as plt



def bs(path='c:\users\sammy\desktop\NanoNet\input_samples', kk=0, flag=True):
    """
    This function computes the band gap / band structure for Silicon Nanowire
    :param path: directory path to where xyz input files are stored
    :param flag: boolean statements
    :return: band gap / band structure
    """
    # define orbitals sets
    tb.Atom.orbital_sets = {'Si': 'SiliconSP3D5S', 'H': 'HydrogenS'}
    band_gaps = []
    band_structures = []
    for xyz_file in os.listdir(path):
        if xyz_file.endswith('xyz'):

            widths = [s for s in xyz_file if s.isdigit()]
            width = int(''.join(widths))
            print(width)

            hamiltonian = tb.Hamiltonian(xyz=os.path.join(path, xyz_file), nn_distance=2.4)
            hamiltonian.initialize()

            if flag:
                plt.axis('off')
                plt.imshow(np.log(np.abs(hamiltonian.h_matrix)))
                plt.savefig('hamiltonian.pdf')
                plt.show()

            a_si = 5.50
            PRIMITIVE_CELL = [[0, 0, a_si]]
            hamiltonian.set_periodic_bc(PRIMITIVE_CELL)

            num_points = len(kk)

            band_structure = []

            for jj in xrange(num_points):
                print(jj)
                vals, _ = hamiltonian.diagonalize_periodic_bc([0, 0, kk[jj]])
                band_structure.append(vals)

            band_structure = np.array(band_structure)
            band_structures.append(band_structure)

            cba = band_structure.copy()
            vba = band_structure.copy()

            cba[cba < 0] = 1000
            vba[vba > 0] = -1000

            band_gap = np.min(cba) - np.max(vba)
            band_gaps.append(band_gap)

            if flag:
                fig1, ax1 = plt.subplots()
                ax1.axes.set_xlim(0.0, 0.5)
                #ax1.axes.set_ylim(0.00, -0.04)
                ax1.axes.set_ylim(-0.3, 0.0)
                ax1.plot(kk, np.sort(np.real(vba)))
                ax1.set_xlabel(r'Wave vector ($\frac{\pi}{a}$)')
                ax1.set_ylabel(r'Energy (eV)')
                ax1.set_title('Valence band, NW width={} u.c.'.format(width))
                fig1.tight_layout()
                plt.savefig('bs_vb.pdf')

                fig2, ax2 = plt.subplots()
                ax1.axes.set_xlim(0.0, 0.5)
                #ax2.axes.set_ylim(0.0, 0.06)
                ax2.axes.set_ylim(2.0, 2.5)
                ax2.plot(kk, np.sort(np.real(cba)))
                ax2.set_xlabel(r'Wave vector ($\frac{\pi}{a}$)')
                ax2.set_ylabel(r'Energy (eV)')
                ax2.set_title('Conduction band, NW width={} u.c.'.format(width))
                fig2.tight_layout()
                plt.savefig('bs_cb.pdf')

    return band_gaps, band_structures
