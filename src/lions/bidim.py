"""2D histograms (bidims) in the GANIL binary format."""

import struct
import sys

import numpy as np

from .constants import C_LIGHT
from .beam import phase_space_coords
from .fortran_text import fortran_float, parse_floats
from .io import prompt
from .state import TrapCircle


class BidimWriter:
    """Writer for up to 50 bidim files, each buffered in blocks of 4096 int32 events."""

    def __init__(self):
        self.n_open = 0
        self.first = 0
        self.last = 0
        self.x_bits = np.zeros(50, dtype=np.int32)
        self.y_bits = np.zeros(50, dtype=np.int32)
        self.x_max = np.zeros(50, dtype=np.int32)
        self.y_max = np.zeros(50, dtype=np.int32)
        self.fill = np.zeros(50, dtype=np.int32)  # next free buffer slot, 1-based as in the Fortran
        self.buffer = np.zeros((4096, 50), dtype=np.int32)
        self.files = {}

    def open(self, count, first, x_bits, y_bits, titles, x_labels, y_labels, bounds):
        """Create files <first>.bid ... <first+count-1>.bid and write their headers. Returns 0 or an error code."""
        self.n_open = 0
        date = [0, 0, 0]
        time = [0, 0, 0]
        # Stored header strings, kept as in the GANIL format
        kind = 'bidim'
        place = 'ganil'
        author = 'inconnu'
        fmt = 'Binaire'
        padding = ' ' * 120

        if count <= 0 or count > 50:
            print('bidim open: wrong number of bidims')
            return 1
        if first <= 0:
            print('bidim open: wrong first bidim number')
            return 2

        self.first = first
        self.last = first + count - 1

        for i in range(count):
            mx = 2**x_bits[i] - 1
            if mx > 8191:
                print(f'bidim open: bidim {i+1}, wrong number of x bits: {x_bits[i]}')
                return 3
            self.x_bits[i] = x_bits[i]
            self.x_max[i] = mx
            my = 2**y_bits[i] - 1
            if my > 8191:
                print(f'bidim open: bidim {i+1}, wrong number of y bits: {y_bits[i]}')
                return 4
            self.y_bits[i] = y_bits[i]
            self.y_max[i] = my

        for i in range(count):
            num = first + i
            path = f'{num}.bid'
            print(f'bidim file: {path}')
            try:
                f = open(path, 'wb')
                self.files[i] = f
                maxx = self.x_max[i]
                maxy = self.y_max[i]
                xmin = f'{bounds[0, 0, i]:.10e}'
                xmax = f'{bounds[1, 0, i]:.10e}'
                ymin = f'{bounds[0, 1, i]:.10e}'
                ymax = f'{bounds[1, 1, i]:.10e}'
                header = struct.pack('i', num)
                header += kind.encode('ascii').ljust(10)
                header += struct.pack('3i', *date)
                header += struct.pack('3i', *time)
                header += place.encode('ascii').ljust(10)
                header += author.encode('ascii').ljust(10)
                header += titles[i].encode('ascii').ljust(25)
                header += x_labels[i].encode('ascii').ljust(25)
                header += y_labels[i].encode('ascii').ljust(25)
                header += struct.pack('2i', maxx, maxy)
                header += fmt.encode('ascii').ljust(10)
                header += xmin.encode('ascii').ljust(20)
                header += xmax.encode('ascii').ljust(20)
                header += ymin.encode('ascii').ljust(20)
                header += ymax.encode('ascii').ljust(20)
                header += padding.encode('ascii')[:120]
                f.write(header)
                self.fill[i] = 1
            except IOError:
                print(f'bidim open: cannot open bidim {i+1}')
                for k in range(i):
                    self.files[k].close()
                return 5

        print('bidim open: all bidim files opened')
        self.n_open = count
        return 0

    def add(self, n, ix, iy):
        """Add one event (ix, iy) to bidim number n. Returns 0 or 1 if n is out of range."""
        if n < self.first or n > self.last:
            print(f'bidim add: wrong bidim number: {n}')
            print(f'...outside: {self.first} {self.last}')
            return 1
        num = n - self.first
        ix = max(0, min(ix, self.x_max[num]))
        iy = max(0, min(iy, self.y_max[num]))
        if self.fill[num] == 4097:
            self._flush(num)
            self.fill[num] = 1
        idx = self.fill[num] - 1
        self.buffer[idx, num] = ix * 65536 + iy
        self.fill[num] += 1
        return 0

    def close(self):
        """Flush the buffers, with a 9999/9999 end marker after a partial block, and close the files."""
        for num in range(self.n_open):
            if self.fill[num] == 4097:
                self._flush(num)
                self.fill[num] = 1
            elif self.fill[num] != 1:
                idx = self.fill[num] - 1
                self.buffer[idx, num] = 9999 * 65536 + 9999
                self._flush(num)
            if num in self.files:
                self.files[num].close()
        return 0

    def _flush(self, num):
        """Write the full 4096-event buffer of bidim num (0-based)."""
        if num in self.files:
            data = self.buffer[:4096, num]
            self.files[num].write(data.tobytes())


def setup_bidims(d, bidim, numbid=None):
    """Set up the bidim diagnostics. numbid=None: ask on stdin."""
    d.bunch.in_window[:] = True
    d.bidims.every = 1
    if numbid is None:
        numbid = int(prompt('first bidim number (0 for no bidim):'))
    d.bidims.first = numbid
    if d.bidims.first == 0:
        return

    d.bidims.every = int(prompt('bidim period in time steps:'))
    d.bidims.centre_mode = prompt('plot centres of curvature or barycentres (c or b):')
    d.bidims.bary_size = abs(d.run.steps_per_turn // d.bidims.every)
    d.bidims.bary_count = 0

    vals = parse_floats(prompt(
        'bidim box (xmin xmax ymin ymax in metres):'))
    box = np.zeros((2, 2))
    box[0, 0] = vals[0]
    box[0, 1] = vals[1]
    box[1, 0] = vals[2]
    box[1, 1] = vals[3]

    x_bits = np.full(20, 13, dtype=np.int32)
    y_bits = np.full(20, 13, dtype=np.int32)
    titles = [''] * 20
    x_labels = [''] * 20
    y_labels = [''] * 20
    bounds = np.zeros((2, 2, 50), dtype=np.float64)

    titles[0] = 'ion trajectory'
    x_labels[0] = 'x'
    y_labels[0] = 'y'
    if d.bidims.centre_mode == 'c':
        titles[1] = 'centres of curvature'
    else:
        titles[1] = 'barycentres'
    x_labels[1] = 'x'
    y_labels[1] = 'y'

    d.bidims.start_step = 1
    if d.bunch.n_particles > 1:
        titles[2] = 'ion trajectories'
        titles[3] = 'correlation dphi-dw/w'
        titles[4] = 'correlation r-dpr/p0'
        vals = parse_floats(prompt('bidim dphi-dw/w: dphi_max dw/w_max:'))
        dphi_max = vals[0]
        dw_max = vals[1]
        vals = parse_floats(prompt('bidim dr-dpr/p0: dr_max dpr/p0_max:'))
        dr_max = vals[0]
        dpr_max = vals[1]
        vals = parse_floats(prompt('bidim dz-dpz/p0: dz_max dpz/p0_max:'))
        dz_max = vals[0]
        dpz_max = vals[1]
        print(f'number of steps: {d.run.n_steps}')
        d.bidims.start_step = int(prompt('first step stored in the bidims (1 to number of steps):'))

    vals = parse_floats(prompt('barycentre window in metres (xmin xmax ymin ymax):'))
    bary_xmin, bary_xmax, bary_ymin, bary_ymax = vals[:4]

    bounds[0, 0, 0] = box[0, 0]
    bounds[1, 0, 0] = box[0, 1]
    bounds[0, 1, 0] = box[1, 0]
    bounds[1, 1, 0] = box[1, 1]
    bounds[0, 0, 1] = bary_xmin
    bounds[1, 0, 1] = bary_xmax
    bounds[0, 1, 1] = bary_ymin
    bounds[1, 1, 1] = bary_ymax

    if d.bunch.n_particles > 1:
        n_bidims = 17
    else:
        n_bidims = 5

    err = bidim.open(n_bidims, d.bidims.first, x_bits, y_bits, titles, x_labels, y_labels, bounds)
    if err != 0:
        print(f'bidim open failed, error {err}')
        sys.exit(1)

    d.bidims.ax1[0] = 8191.0 / (box[0, 1] - box[0, 0])
    d.bidims.ax0[0] = -d.bidims.ax1[0] * box[0, 0]
    d.bidims.ay1[0] = 8191.0 / (box[1, 1] - box[1, 0])
    d.bidims.ay0[0] = -d.bidims.ay1[0] * box[1, 0]
    d.bidims.ax1[1] = 8191.0 / (bary_xmax - bary_xmin)
    d.bidims.ax0[1] = -d.bidims.ax1[1] * bary_xmin
    d.bidims.ay1[1] = 8191.0 / (bary_ymax - bary_ymin)
    d.bidims.ay0[1] = -d.bidims.ay1[1] * bary_ymin

    if d.bunch.n_particles > 1:
        d.bidims.ax1[3] = 8191.0 / 2.0 / dphi_max
        d.bidims.ax0[3] = d.bidims.ax1[3] * dphi_max
        d.bidims.ay1[3] = 8191.0 / 2.0 / dw_max
        d.bidims.ay0[3] = d.bidims.ay1[3] * dw_max
        d.bidims.ax1[4] = 8191.0 / 2.0 / dr_max
        d.bidims.ax0[4] = d.bidims.ax1[4] * dr_max
        d.bidims.ay1[4] = 8191.0 / 2.0 / dpr_max
        d.bidims.ay0[4] = d.bidims.ay1[4] * dpr_max
        d.bidims.ax1[5] = 8191.0 / 2.0 / dr_max
        d.bidims.ax0[5] = d.bidims.ax1[5] * dr_max
        d.bidims.ay1[5] = 8191.0 / 2.0 / dphi_max
        d.bidims.ay0[5] = d.bidims.ay1[5] * dphi_max
        d.bidims.ax1[6] = 8191.0 * d.bidims.every / d.run.n_steps
        d.bidims.ax0[6] = 0.0
        d.bidims.ay1[6] = 8191.0 / 2.0 / dr_max
        d.bidims.ay0[6] = d.bidims.ay1[6] * dr_max
        d.bidims.ax1[7] = 8191.0 / 2.0 / dz_max
        d.bidims.ax0[7] = d.bidims.ax1[7] * dz_max
        d.bidims.ay1[7] = 8191.0 / 2.0 / dpz_max
        d.bidims.ay0[7] = d.bidims.ay1[7] * dpz_max
        d.bidims.ax1[8] = 8191.0 * d.bidims.every / d.run.n_steps
        d.bidims.ax0[8] = 0.0
        d.bidims.ay1[8] = 8191.0 / 2.0 / dz_max
        d.bidims.ay0[8] = d.bidims.ay1[8] * dz_max

    # Marker file: lines and circles (pillars, traps) drawn into bidim 1, in mm.
    # Circles with flag < -10 are ion traps; flag <= -9000 means trap only, not drawn.
    marker_file = prompt('marker file (lines, pillars, etc.), or no:')
    d.traps.circles = []
    if marker_file != 'no':
        try:
            with open(marker_file, 'r') as f:
                for _ in range(3):  # display window, grid size, niter omer: unused
                    f.readline()
                n_points = 360
                while True:
                    kind = f.readline().strip()
                    if kind == 'fintyp':
                        break
                    elif kind == 'droite':
                        n_lines = int(f.readline().split()[0])
                        for _ in range(n_lines):
                            vals = f.readline().split()
                            x1, y1, x2, y2 = [float(v) for v in vals[:4]]
                            step_x = (x2 - x1) / (n_points - 1.0)
                            step_y = (y2 - y1) / (n_points - 1.0)
                            for kk in range(n_points):
                                xd = x1 + kk * step_x
                                yd = y1 + kk * step_y
                                ix = int(d.bidims.ax1[0] * xd / 1000.0 + d.bidims.ax0[0])
                                iy = int(d.bidims.ay1[0] * yd / 1000.0 + d.bidims.ay0[0])
                                bidim.add(d.bidims.first, ix, iy)
                    elif kind == 'cercle':
                        n_circles = int(f.readline().split()[0])
                        for _ in range(n_circles):
                            vals = f.readline().split()
                            xc, yc, rc = fortran_float(vals[0]), fortran_float(vals[1]), fortran_float(vals[2])
                            theta1 = fortran_float(vals[3]) * np.pi / 180.0
                            theta2 = fortran_float(vals[4]) * np.pi / 180.0
                            flag = fortran_float(vals[6])
                            dtheta = (theta2 - theta1) / (n_points - 1.0)
                            if flag < -10.0:
                                d.traps.circles.append(TrapCircle(xc / 1000.0, yc / 1000.0, rc / 1000.0))
                            if flag > -9000.0:
                                for kk in range(n_points):
                                    theta = theta1 + kk * dtheta
                                    xd = xc + rc * np.cos(theta)
                                    yd = yc + rc * np.sin(theta)
                                    ix = int(d.bidims.ax1[0] * xd / 1000.0 + d.bidims.ax0[0])
                                    iy = int(d.bidims.ay1[0] * yd / 1000.0 + d.bidims.ay0[0])
                                    bidim.add(d.bidims.first, ix, iy)
                    else:
                        print(f'marker file: unknown keyword: {kind}')
                        sys.exit(1)
        except Exception as e:
            print(f'error reading the marker file: {e}')

    if len(d.traps.circles) != 0:
        print('warning: ion traps are defined')
        print(f'number of traps: {len(d.traps.circles)}')
        d.traps.last_step = int(prompt('last step with active traps:'))
        print('trapped particles: (0 or 1)')
        print('0   : stay stuck in the traps')
        print('1   : keep going, merged with the central particle')
        d.traps.follow = int(prompt(''))


def fill_bidims(d, bidim):
    """Add the current step to the bidims."""
    vec1 = np.zeros(d.bunch.n_particles)
    vec2 = np.zeros(d.bunch.n_particles)
    vec3 = np.zeros(d.bunch.n_particles)
    vec4 = np.zeros(d.bunch.n_particles)
    vec5 = np.zeros(d.bunch.n_particles)
    vec6 = np.zeros(d.bunch.n_particles)
    if d.run.step < d.bidims.start_step:
        return

    for k in range(d.bunch.n_particles):
        if not d.bunch.in_window[k]:
            continue
        ix = int(d.bidims.ax1[0] * d.bunch.x[k] + d.bidims.ax0[0])
        iy = int(d.bidims.ay1[0] * d.bunch.y[k] + d.bidims.ay0[0])
        if 0 <= ix <= 8191 and 0 <= iy <= 8191:
            if k == 0:
                bidim.add(d.bidims.first, ix, iy)
                if d.bunch.n_particles > 1:
                    bidim.add(d.bidims.first + 2, ix, iy)
            else:
                bidim.add(d.bidims.first + 2, ix, iy)
        if k == 0:
            xc = d.bunch.x[k] + d.bunch.py[k] * d.run.dt / d.bunch.boris_e_coef[k] / d.bunch.bz[k]
            yc = d.bunch.y[k] - d.bunch.px[k] * d.run.dt / d.bunch.boris_e_coef[k] / d.bunch.bz[k]
            if d.bidims.centre_mode == 'c':
                ix2 = int(d.bidims.ax1[0] * xc + d.bidims.ax0[0])
                iy2 = int(d.bidims.ay1[0] * yc + d.bidims.ay0[0])
                if 0 <= ix2 <= 8191 and 0 <= iy2 <= 8191:
                    bidim.add(d.bidims.first + 1, ix2, iy2)
            if d.bidims.centre_mode == 'b':
                d.bidims.bary_count += 1
                d.bidims.bary_x[d.bidims.bary_count - 1] = xc
                d.bidims.bary_y[d.bidims.bary_count - 1] = yc
                if d.bidims.bary_count == d.bidims.bary_size:
                    d.bidims.bary_count = 0
                    xbary = np.sum(d.bidims.bary_x[:d.bidims.bary_size]) / d.bidims.bary_size
                    ybary = np.sum(d.bidims.bary_y[:d.bidims.bary_size]) / d.bidims.bary_size
                    ix2 = int(d.bidims.ax1[1] * xbary + d.bidims.ax0[1])
                    iy2 = int(d.bidims.ay1[1] * ybary + d.bidims.ay0[1])
                    if 0 <= ix2 <= 8191 and 0 <= iy2 <= 8191:
                        bidim.add(d.bidims.first + 1, ix2, iy2)

    if d.bunch.n_particles > 1 and len(d.bunch.families) == 1:
        phase_space_coords(d.bunch)
        k = 0
        nbr = d.bunch.families[k].size
        i1 = d.bunch.families[k].start
        gamma_ref = np.sqrt(1.0 + d.bunch.p[i1]**2)
        vec1[:nbr] = d.bunch.ps[4, :nbr] * d.rf.omega / C_LIGHT / d.bunch.p[i1] * 180.0 / np.pi * gamma_ref
        vec2[:nbr] = d.bunch.ps[5, :nbr] * d.bunch.p[i1] / gamma_ref / (gamma_ref - 1.0)
        vec3[:nbr] = d.bunch.ps[0, :nbr] * 1000.0
        vec4[:nbr] = d.bunch.ps[1, :nbr] / d.bunch.p[i1] * 1000.0
        vec5[:nbr] = d.bunch.ps[2, :nbr] * 1000.0
        vec6[:nbr] = d.bunch.ps[3, :nbr] / d.bunch.p[i1] * 1000.0

        for i in range(nbr):
            if not d.bunch.in_window[i]:
                continue
            ix = int(d.bidims.ax1[3] * vec1[i] + d.bidims.ax0[3])
            iy = int(d.bidims.ay1[3] * vec2[i] + d.bidims.ay0[3])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.add(d.bidims.first + 3, ix, iy)

            ix = int(d.bidims.ax1[4] * vec3[i] + d.bidims.ax0[4])
            iy = int(d.bidims.ay1[4] * vec4[i] + d.bidims.ay0[4])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.add(d.bidims.first + 4, ix, iy)

            ix = int(d.bidims.ax1[7] * vec5[i] + d.bidims.ax0[7])
            iy = int(d.bidims.ay1[7] * vec6[i] + d.bidims.ay0[7])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.add(d.bidims.first + 7, ix, iy)

        for i in range(nbr):
            if not d.bunch.in_window[i]:
                continue
            ix = int(d.bidims.ax1[5] * vec3[i] + d.bidims.ax0[5])
            iy = int(d.bidims.ay1[5] * vec1[i] + d.bidims.ay0[5])
            ix = max(0, min(ix, 8191))
            iy = max(0, min(iy, 8191))
            bidim.add(d.bidims.first + 5, ix, iy)

        ix = int(d.bidims.ax1[6] * d.run.step / d.bidims.every + d.bidims.ax0[6])
        ix = max(0, min(ix, 8191))
        for i in range(nbr):
            if not d.bunch.in_window[i]:
                continue
            iy = int(d.bidims.ay1[6] * vec3[i] + d.bidims.ay0[6])
            iy = max(0, min(iy, 8191))
            bidim.add(d.bidims.first + 6, ix, iy)
            iy = int(d.bidims.ay1[8] * vec5[i] + d.bidims.ay0[8])
            iy = max(0, min(iy, 8191))
            bidim.add(d.bidims.first + 8, ix, iy)
