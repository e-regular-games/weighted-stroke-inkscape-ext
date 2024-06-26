#!/usr/bin/env python
# coding=utf-8
#
# Copyright (C) 2024, E-Regular Games, e.regular.games.llc@gmail.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
Description of this extension
"""

import math
import cmath
import numpy as np

import inkex
from inkex.paths import (
    Move,
    Line,
    Curve,
    ZoneClose,
    Arc,
    Path,
    Vert,
    Horz,
    TepidQuadratic,
    Quadratic,
    Smooth,
)
from inkex.transforms import Vector2d
from inkex.bezier import beziertatslope, beziersplitatt, bezierlength, bezierslopeatt, bezierparameterize, bezierpointatt, maxdist

X=0
Y=1

M=0 # y = Mx+b
B=1 # y = Mx+b

TOL = 0.001 # same precision shown in the inkscape ui

RES = 4 # pair with the Q below

# rats(inv(M*T)) where T = [t=0, t=1/3, t=2/3, t=1]
# note P*Q = C, where C is the control points of the bezier
# and P are the points on the curve. 
Q = np.array([[1.0, -5.0/6.0,  1.0/3.0, 0.0],
              [0.0,      3.0, -3.0/2.0, 0.0],
              [0.0, -3.0/2.0,      3.0, 0.0],
              [0.0,  1.0/3.0, -5.0/6.0, 1.0]])

def root_wrapper_3(root_a, root_b, root_c, root_d, dbg):
    roots = np.zeros(root_a.shape + (3,), dtype=np.cdouble)
    
    # Monics formula, see
    # http://en.wikipedia.org/wiki/Cubic_function#Monic_formula_of_roots
    mono_a, mono_b, mono_c = (root_b / root_a, root_c / root_a, root_d / root_a)
    m = 2.0 * np.power(mono_a,3) - 9.0 * mono_a * mono_b + 27.0 * mono_c
    k = np.power(mono_a,2) - 3.0 * mono_b
    n = np.power(m,2) - 4.0 * np.power(k,3)
    w1 = -0.5 + 0.5 * np.emath.sqrt(-3.0)
    w2 = -0.5 - 0.5 * np.emath.sqrt(-3.0)
    
    m1 = np.zeros(m.shape, dtype=np.cdouble)
    n1 = np.zeros(m.shape, dtype=np.cdouble)
    
    # n < 0 cases
    idx_n_imag = n < 0
    m1[idx_n_imag] = np.power((m[idx_n_imag] + np.emath.sqrt(n[idx_n_imag])) / 2, 1.0 / 3.0)
    n1[idx_n_imag] = np.power((m[idx_n_imag] - np.emath.sqrt(n[idx_n_imag])) / 2, 1.0 / 3.0)
    
    # n > 0 cases
    idx_n_real = np.logical_not(idx_n_imag)
    sel = np.full(m.shape, False)
    sel[idx_n_real] = (m[idx_n_real] + np.sqrt(n[idx_n_real])) < 0
    m1[sel] = -1 * np.power(-1 * (m[sel] + np.sqrt(n[sel])) / 2, 1.0 / 3)
    sel = np.full(m.shape, False)
    sel[idx_n_real] = (m[idx_n_real] + np.sqrt(n[idx_n_real])) >= 0
    m1[sel] = np.power((m[sel] + np.sqrt(n[sel])) / 2, 1.0 / 3)
    sel = np.full(m.shape, False)
    sel[idx_n_real] = (m[idx_n_real] - np.sqrt(n[idx_n_real])) < 0
    n1[sel] = -1 * np.power(-1 * (m[sel] - np.sqrt(n[sel])) / 2, 1.0 / 3)
    sel = np.full(m.shape, False)
    sel[idx_n_real] = (m[idx_n_real] - np.sqrt(n[idx_n_real])) >= 0
    n1[sel] = np.power((m[sel] - np.sqrt(n[sel])) / 2, 1.0 / 3)

    roots[...,0] = -1.0 / 3 * (mono_a + m1 + n1)
    roots[...,1] = -1.0 / 3 * (mono_a + w1 * m1 + w2 * n1)
    roots[...,2] = -1.0 / 3 * (mono_a + w2 * m1 + w1 * n1)

    valid = np.logical_and(np.imag(roots) == 0, np.logical_and(np.real(roots) >= 0, np.real(roots) <= 1))
    roots[np.logical_not(valid)] = 0.0
    
    return np.real(roots), valid

def root_wrapper_2(root_b, root_c, root_d):
    roots = np.zeros(root_b.shape + (2,), dtype=np.cdouble)
    valid = np.full(roots.shape, False)
    det = np.power(root_c,2) - 4.0 * root_b * root_d

    roots_2 = det != 0
    roots[roots_2,0] = (-root_c[roots_2] + np.emath.sqrt(det[roots_2])) / (2.0 * root_b[roots_2])
    roots[roots_2,1] = (-root_c[roots_2] - np.emath.sqrt(det[roots_2])) / (2.0 * root_b[roots_2])
    valid[roots_2,0:2] = 1
    
    roots_1 = det == 0
    roots[roots_1,0] = -root_c[roots_1] / (2.0 * root_b[roots_1])
    valid[roots_1,0:1] = 1

    valid = np.logical_and(valid, np.logical_and(np.imag(roots) == 0, np.logical_and(np.real(roots) >= 0, np.real(roots) <= 1)))
    roots[np.logical_not(valid)] = 0.0
    
    return roots, valid

def root_wrapper_1(root_c, root_d):
    roots = -root_d / root_c
    valid = np.logical_and(roots >= 0, roots <= 1)
    roots[np.logical_not(valid)] = 0.0
    return roots.reshape((-1,1)), valid.reshape((-1,1))

def root_wrapper(root_a, root_b, root_c, root_d, dbg):
    """Get the Cubic function, moic formular of roots, simple root"""

    num_roots = np.zeros(root_a.shape)
    
    num_roots[root_a != 0] = 3
    num_roots[np.logical_and(num_roots == 0, root_b != 0)] = 2
    num_roots[np.logical_and(num_roots == 0, root_c != 0)] = 1

    roots = np.zeros(root_a.shape + (3,))
    valid = np.full(root_a.shape + (3,), False)

    roots_3 = num_roots == 3
    if np.any(roots_3):
        roots[roots_3,0:3], valid[roots_3,0:3] = \
            root_wrapper_3(root_a[roots_3], root_b[roots_3], root_c[roots_3], root_d[roots_3], dbg)

    roots_2 = num_roots == 2
    if np.any(roots_2):
        roots[roots_2,0:2], valid[roots_2,0:2] = \
            root_wrapper_2(root_b[roots_2], root_c[roots_2], root_d[roots_2])
    
    roots_1 = num_roots == 1
    if np.any(roots_1):
        roots[roots_1,0:1], valid[roots_1,0:1] = \
            root_wrapper_1(root_c[roots_1], root_d[roots_1])

    return roots, valid

class WeightedStroke(inkex.EffectExtension):
    """Please rename this class, don't keep it unnamed"""
    def add_arguments(self, pars):
        pars.add_argument("--weight", type=float,\
            help="The boldness percentage to increase or decrease the path weight.")
        pars.add_argument("--debug", type=inkex.Boolean, help="Render debug lines and points.", default=False)
        pars.add_argument("--shrink", type=inkex.Boolean, help="Shrink the shape.", default=False)

        
    def _to_cubics(self, beziers):
        cubics = np.zeros(beziers.shape)
        cubics[:,:,3] = beziers[:,:,0]
        cubics[:,:,2] = 3 * (beziers[:,:,1] - beziers[:,:,0])
        cubics[:,:,1] = 3 * (beziers[:,:,2] - beziers[:,:,1]) - cubics[:,:,2]
        cubics[:,:,0] = beziers[:,:,3] - beziers[:,:,0] - cubics[:,:,2] - cubics[:,:,1]
        return cubics

    def dist(self, x0, y0, x1, y1):
        return math.sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1))

    # turn all paths into bezier curves, including lines, verticals, horizontals, and half-line/half-curve.
    def _curves_matrix(self, abs_path):
        reset_origin = True
        beziers = np.zeros((2, len(abs_path), 4))
        linear = np.full((len(abs_path),), False)
        i = 0
        for cmd_proxy in abs_path.proxy_iterator():
            # for each subpath, reset the origin of the following computations to the first
            # node of the subpath -> i.e. after a Z command, move the origin to the end point
            # of the next command
            if reset_origin:
                first_point = cmd_proxy.end_point
                reset_origin = False
            if isinstance(cmd_proxy.command, ZoneClose):
                #todo: handle case where ZoneClose is a striaght line from end of bezier to starting point
                reset_origin = True
            if isinstance(cmd_proxy.command, (Curve, Smooth, TepidQuadratic, Quadratic, Arc)):
                prev = cmd_proxy.previous_end_point
                for curve in cmd_proxy.to_curves():
                    bez = curve.to_bez()
                    d = self.dist(bez[2][0], bez[2][1], prev.x, prev.y)
                    if d > TOL:
                        beziers[X,i] = [prev.x, bez[0][0], bez[1][0], bez[2][0]]
                        beziers[Y,i] = [prev.y, bez[0][1], bez[1][1], bez[2][1]]
                        i += 1
            elif isinstance(cmd_proxy.command, (Line, ZoneClose, Vert, Horz)):
                prev = cmd_proxy.previous_end_point
                n = None
                if isinstance(cmd_proxy.command, (ZoneClose,)):
                    n = (first_point.x, first_point.y)
                elif isinstance(cmd_proxy.command, (Vert,)):
                    n = (prev.x, cmd_proxy.command.args[0])
                elif isinstance(cmd_proxy.command, (Horz,)):
                    n = (cmd_proxy.command.args[1], prev.y)
                else:
                    n = cmd_proxy.command.args
                if self.dist(prev.x, prev.y, n[0], n[1]) > TOL:
                    beziers[X,i] = [
                        prev.x,
                        2.0*prev.x/3.0 + n[0]/3.0,
                        prev.x/3.0 + 2.0*n[0]/3.0,
                        n[0]
                    ]
                    beziers[Y,i] = [
                        prev.y,
                        2.0*prev.y/3.0 + n[1]/3.0,
                        prev.y/3.0 + 2.0*n[1]/3.0,
                        n[1]
                    ]
                    linear[i] = True
                    i += 1
        return beziers[:,:i], linear[:i]

    def _slope_at_t(self, beziers, cubics, t):
        # generic result, good for t != 0 and t != 1.
        m = 3 * cubics[:,:,0:1] * np.power(t, 2) + 2 * cubics[:,:,1:2] * t + cubics[:,:,2:3]
        
        replace_m0 = t == 0
        if np.any(replace_m0):
            pf = beziers[:,:,1:2]
            fixed = np.all(beziers[:,:,1:2] == beziers[:,:,0:1], axis=0)[:,0]
            pf[:,fixed] = beziers[:,fixed,2:3]
            m[:,:,replace_m0] = pf - beziers[:,:,0:1]

            replace_m1 = t == 1
        if np.any(replace_m1):
            pi = beziers[:,:,2:3]
            fixed = np.all(beziers[:,:,3:4] == beziers[:,:,2:3], axis=0)[:,0]
            pi[:,fixed] = beziers[:,fixed,1:2]
            m[:,:,replace_m1] = beziers[:,:,3:4] - pi

        return m
        
    def _point_at_t(self, cubics, t):
        return cubics[:,:,0:1] * np.power(t, 3) \
            + cubics[:,:,1:2] * np.power(t, 2) \
            + cubics[:,:,2:3] * t \
            + cubics[:,:,3:4]

    def _slope_coefficients(self, slopes):
        # slopes -> [x/y, N]
        def_slope = slopes[Y,:] != 0
        inf_slope = np.logical_and(np.logical_not(def_slope), slopes[X,:] != 0)
        dne = np.logical_and(slopes[X,:] == 0, slopes[Y,:] == 0)
        coef1 = np.ones((1,slopes.shape[1],))
        coef2 = np.ones((1,slopes.shape[1],))
        
        coef1[0,def_slope] = slopes[X,def_slope] / slopes[Y,def_slope]
        coef1[0,dne] = 0.0
        coef2[0,inf_slope] = slopes[Y,inf_slope] / slopes[X,inf_slope]
        coef2[0,dne] = 0.0

        return coef1, coef2
        
    
    def _intersections_closest(self, lines, cubics):
        """
        Where each line (n) intersects with each cubic (m).
        uses: numpy arrays for optimized calculations.
        lines is:
          [
            [[slope0_x, point0_x], [slope1_x, point1_x], ...],
            [[slope0_y, point0_y], [slope1_y, point1_y], ...],
          ]
        
        """

        I = 3 # intersects
        coef1, coef2 = self._slope_coefficients(lines[...,M])
        
        # cubic intersection coefficients
        a = coef1 * cubics[Y,:,0:1] - coef2 * cubics[X,:,0:1]
        b = coef1 * cubics[Y,:,1:2] - coef2 * cubics[X,:,1:2]
        c = coef1 * cubics[Y,:,2:3] - coef2 * cubics[X,:,2:3]
        d = coef1 * (cubics[Y,:,3:4] - lines[Y,:,B:B+1].transpose()) \
            - coef2 * (cubics[X,:,3:4] - lines[X,:,B:B+1].transpose())

        # of the shape (cubic(m), line(n), intersect(3)]
        roots, valid = root_wrapper(a, b, c, d, self)

        # determine the distance from each point to each intersect.
        # use argmin to determine the shortest distance

        on_line_exp = (np.ones((1,I)) * np.arange(lines.shape[1]).reshape((lines.shape[1],1))).reshape((1,lines.shape[1]*I))
        on_line_exp = (on_line_exp * np.ones((cubics.shape[1],1))).reshape((cubics.shape[1]*lines.shape[1]*I,)).astype(np.int32)
        
        int_cubic = (np.ones((1,roots.shape[1] * I)) * np.arange(cubics.shape[1]).reshape((cubics.shape[1],1))) \
            .reshape((cubics.shape[1] * roots.shape[1] * I,)).astype(np.int32)

        sort_idx = np.argsort(on_line_exp, axis=0)
        on_line_exp = on_line_exp[sort_idx]
        int_cubic = int_cubic[sort_idx]

        roots_exp = roots.reshape((roots.shape[0] * roots.shape[1] * I,))[sort_idx]
        valid_exp = valid.reshape((valid.shape[0] * valid.shape[1] * I,))[sort_idx]
        cubics_exp = cubics[:,int_cubic]

        intersects = np.zeros((2, roots_exp.shape[0]))
        intersects[:,valid_exp] = self._point_at_t(cubics_exp[:,valid_exp], roots_exp[valid_exp].reshape(-1,1))[:,:,0]
        points = lines[:,on_line_exp,B]

        dist = np.zeros((roots_exp.shape[0],))
        diff = intersects[:,valid_exp] - points[:,valid_exp]

        lines_slope_sign = np.sign(lines[:,:,M])
        diff_slope_sign = np.sign(diff)
        inside_shape = np.full(valid_exp.shape, False)
        inside_shape[valid_exp] = \
            np.logical_and(diff_slope_sign[0] == lines_slope_sign[0,on_line_exp[valid_exp]],
                           diff_slope_sign[1] == lines_slope_sign[1,on_line_exp[valid_exp]])

        dist[valid_exp] = np.sqrt(diff[0]*diff[0] + diff[1]*diff[1])
        valid_exp = np.logical_and(np.logical_and(valid_exp, dist > TOL), inside_shape)

        valid_sort_shaped = (1.0 - valid_exp).reshape((lines.shape[1], I * cubics.shape[1]))
        dist_shaped = dist.reshape((lines.shape[1], I * cubics.shape[1]))
        closest_index = I * cubics.shape[1] * np.arange(0, lines.shape[1]) + np.lexsort((dist_shaped, valid_sort_shaped))[:,0]

        closest = np.full(valid_exp.shape, False)
        closest[closest_index] = True
        valid_exp = np.logical_and(valid_exp, closest)
        
        return points[:,valid_exp], intersects[:,valid_exp], on_line_exp[valid_exp]

    def _make_line(self, a, b):
        elem = inkex.PathElement()
        elem.path = inkex.Path(
            [
                Move(a[0], a[1]),
                Line(b[0], b[1]),
                ZoneClose(),
            ]
        )
        elem.style.set_color('#FF0000', 'stroke')
        return elem

    def _msg(self, text):
        if self.options.debug:
            self.msg(text)

    # returns an array of angles with [0] being the angle between curves [0] and [1].
    #         angle [n] is the angle between curves [n] and [0]
    def _angle_between_curves(self, cubics):
        a = np.zeros((cubics.shape[1],1))

        dt0 = self._slope_at_t(cubics, 0)[:,:,0].transpose()
        dt1 = self._slope_at_t(cubics, 1)[:,:,0].transpose()

        dt0_u = dt0 / np.linalg.norm(dt0, axis=1).reshape((-1,1))
        dt1_u = dt1 / np.linalg.norm(dt1, axis=1).reshape((-1,1))

        a[:-1] = np.arccos(np.clip(np.sum(dt0_u[1:] * dt1_u[:-1], axis=1), -1.0, 1.0)).reshape((-1, 1))
        a[-1:] = np.arccos(np.clip(np.sum(dt1_u[-1:] * dt0_u[:1], axis=1), -1.0, 1.0)).reshape((-1, 1))

        return math.pi - a

    def _norm(self, v):
        # v is of form [x/y,...]
        mag = np.sqrt(v[X]*v[X] + v[Y]*v[Y])
        v[X] = v[X] / mag
        v[Y] = v[Y] / mag
        return v

    def _generate_ref_lines(self, beziers, cubics, linear):
        t = np.linspace(0.0, 1.0, num=RES)

        slopes = self._slope_at_t(beziers, cubics, t)
        points = self._point_at_t(cubics, t)

        orthogonals = np.zeros(slopes.shape)
        orthogonals[X] = -1 * slopes[Y]
        orthogonals[Y] = slopes[X]
        orthogonals = self._norm(orthogonals)

        # position 0 refers to the angle between curves 0 and 1.
        flats = np.full((orthogonals.shape[1],1), False)
        flats[:-1,0] = np.all(np.abs(orthogonals[:,:-1,-1] - orthogonals[:,1:,0]) < TOL, axis=0)
        flats[-1,0] = np.all(np.abs(orthogonals[:,-1,-1] - orthogonals[:,0,0]) < TOL, axis=0)
        
        to_calc = np.full((points.shape[1],RES), True)
        to_calc[:, 0] = False
        to_calc[0, 0] = True

        orthogonals[:,:-1,-1] += orthogonals[:,1:,0]
        orthogonals[:,-1,-1] += orthogonals[:,0,0]
        orthogonals[:,0,0] = orthogonals[:,-1,-1]
        orthogonals = self._norm(orthogonals)

        f = np.full((RES,),False)
        f[0] = True
        f[RES-1] = True
        to_calc[linear] = np.logical_and(to_calc[linear], f)
        
        points = points.reshape((2, points.shape[1] * RES))
        on_cubic = (np.ones((1,RES)) * np.arange(0, cubics.shape[1]).reshape((-1,1))).reshape((RES*cubics.shape[1],))
        orthogonals = orthogonals.reshape((2, orthogonals.shape[1] * RES))
        to_calc = to_calc.reshape((to_calc.shape[0] * RES,))

        return points[:,to_calc], orthogonals[:,to_calc], on_cubic[to_calc], flats

    def _outline(self, points, intersects):
        diff = intersects - points
        dist = np.sqrt(diff[0]*diff[0] + diff[1]*diff[1])
        
        delta = -1.0 * diff * self.options.weight / np.max(dist)
        if self.options.shrink:
            delta *= -1
        return points + delta
        
    def effect(self):
        for node in self.svg.selection.filter_nonzero(inkex.PathElement):
            abs_path = node.path.to_absolute()
            beziers, linear = self._curves_matrix(abs_path)
            cubics = self._to_cubics(beziers)

            points, orthogonals, on_cubic, flats = self._generate_ref_lines(beziers, cubics, linear)
            
            lines = np.zeros((2, orthogonals.shape[1], 2)) # [x/y, n, m/b]
            lines[X,:,M] = orthogonals[X]
            lines[Y,:,M] = orthogonals[Y]
            lines[X,:,B] = points[X]
            lines[Y,:,B] = points[Y]

            points, intersects, line_map = self._intersections_closest(lines, cubics)
            
            group = node.getparent().add(inkex.Group())
            points_t = points.transpose()
            intersects_t = intersects.transpose()

            if self.options.debug:
                for p in range(points_t.shape[0]):
                    group.add(self._make_line(points_t[p], intersects_t[p]))

            outline = self._outline(points, intersects)
            elem = inkex.PathElement()
            p = []
            prev = None

            path = self._make_path(cubics.shape[1], on_cubic, outline, flats)
            elem.path = inkex.Path(path)
            elem.style.set_color('#FF0000', 'fill')
            node.getparent().add(elem)

    def _match_slopes(self, s, p):
        wi = np.array([[20, 7], [7, 20]])
        wj = np.linalg.inv([[12, 6], [-6, -12]])
        pi = np.concatenate([p[:,1:2], p[:,2:3]], axis=1)
        pj = np.concatenate([p[:,0:1], p[:,3:4]], axis=1)
        K = np.linalg.inv(s) @ (27 * pi - pj @ wi) @ wj
        return np.concatenate([p[:,0:1], p[:,0:1] + K[0,0] * s[:,0:1], p[:,3:4] - K[1,1] * s[:,1:2], p[:,3:4]], axis=1)

    def _slope(self, prev, points):
        if points.shape[1] == 1:
            s = points[:,0:1] - prev
            return np.concatenate([s, s], axis=1)
        if points.shape[1] == 3:
            # find control points.
            C = np.concatenate([prev, points], axis=1) @ Q
            return np.concatenate([C[:,1:2] - C[:,0:1], C[:,3:4] - C[:,2:3]], axis=1)
        
    def _make_path(self, n, on_cubic, outline, flats):        
        prev_p = None
        p = []

        slopes = np.zeros((n,2,2))
        for i in range(n):
            otl = outline[:,on_cubic == i]
                
            if prev_p is None:
                prev_p = otl[:,0:1]
                otl = otl[:,1:]

            slopes[i,:,:] = self._slope(prev_p, otl)
            prev_p = otl[:,-1:]

        prev_p = None
        for i in range(n):
            otl = outline[:,on_cubic == i]
                
            if prev_p is None:
                prev_p = otl[:,0:1]
                p.append(Move(otl[0,0], otl[1,0]))
                otl = otl[:,1:]

            if otl.shape[1] == 1:
                p.append(Line(otl[0,0], otl[1,0]))
            elif otl.shape[1] == 3:
                # use prev and the 3 points below to compute the bezier curve.
                # the points below are on the curve, not the actual control points.
                s = slopes[i,:,:].copy()
                if flats[(i-1)%n]:
                    s[:,0:1] = (s[:,0:1] + slopes[(i-1)%n,:,1:2]) / 2
                if flats[i]:
                    s[:,1:2] = (s[:,1:2] + slopes[(i+1)%n,:,0:1]) / 2

                P = self._match_slopes(s, np.concatenate([prev_p, otl], axis=1))
                p.append(Curve(P[0,1], P[1,1], P[0,2], P[1,2], P[0,3], P[1,3]))
            prev_p = otl[:,-1:]

        p.append(ZoneClose())
        return p
        
            
if __name__ == '__main__':
    WeightedStroke().run()
