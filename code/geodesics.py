from matplotlib import patches
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import numpy as np
import gmpy2
import math
import os

os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
# use latex (must add to path, found with 'which latex' on terminal and not full path since need other components too, so just texbin)
plt.rcParams.update({
"font.family": "serif",  # use serif/main font for text elements
"text.usetex": True,     # use inline math for ticks
"pgf.rcfonts": False,     # don't setup fonts from rc parameters
"font.serif": "serif",
"text.latex.preamble": r"\usepackage{amssymb, amsmath}"
})

ceil = math.ceil
floor = math.floor
sqrt = math.sqrt

left = -1
bottom = 0
right = 1

def sign(n):
    if n == 0:
        raise ValueError('n cannot be 0')
    if n > 0:
        return 1
    else:
        return -1


def is_square_free(n): 
    if n % 2 == 0: 
        n = n / 2
    if n % 2 == 0: 
        return False

    for i in range(3, int(sqrt(n) + 1)): 
        if n % i == 0: 
            n = n / i 
        if n % i == 0: 
            return False
    return True


# Extract complete convergent (p + sqrt(D))/(2q) of the
# reduced irrational (b+sqrt(D))/(2c), for c > 0, D > 0,
# 0 < b < sqrt(D), sqrt(D) - b < 2c < sqrt(D) + b, 4c | D - b^2
def next_complete_convergent(eta):
    (b, c, D) = eta
    
    if D <= 0:
        raise ValueError('D must be positive', D)
    if (D % 4 != 0) and (D % 4 != 1):
        raise ValueError('D must be 0 or 1 mod 4', D)
    if gmpy2.is_square(D):
        raise ValueError('D cannot be a square', D)

    R = sqrt(D)
    Rfloor = floor(R)

    if c <= 0:
        raise ValueError('c must be positive', c)
    if (b <= 0) or (b > Rfloor):
        raise ValueError('must have 0 < b < R', b, R)
    if (2*c <= Rfloor - b) or (2*c > Rfloor + b):
        raise ValueError('must have R - b < 2c < R + b', R-b, c, R+b)
    if ((D - b*b) % (4 * c) != 0):
        raise ValueError('must have 4c | D - b^2', 4*c, D-b**2)

    x = floor((b + Rfloor)/(2*c))
    p = 2*c*x - b
    a = (b*b - D) // (4*c)
    q = -a + b*x - c*x*x
    return ((p, q, D), x)


# Get cycle of equivalent reduced forms to the given form (a, b, c), and also return
# digits of continued fraction. If cycle of digits is odd then Pell's equation has
# a solution (i.e. there is a unit of norm -1). In this case (a, b, c) ~ (-a, b, -c)
# and must repeat procedure to get full cycle of equivalent forms from continued fraction.
def equiv_reduced_forms(form):
    (a, b, c) = form
    sign_a = sign(a)
    D = b*b - 4*a*c
    cycle_forms = []
    cont_frac_digits = []

    eta = (b, abs(c), D)
    f = eta
    parity = 1
    Pell_has_sol = False
    while True:
        (p, q, DDD) = f
        b = p
        c = (- sign_a) * parity * q
        if (b*b - D) % (4*c) != 0:
            raise ValueError('not a valid form, something went wrong')
        a = (b*b - D) // (4*c)
        cycle_forms.append((a, b, c))

        (f, digit) = next_complete_convergent(f)
        if not Pell_has_sol:
            cont_frac_digits.append(digit)
        parity *= -1
        if f == eta:
            if (parity == -1) and (not Pell_has_sol):
                Pell_has_sol = True
                continue
            else:
                break

    return (cycle_forms, cont_frac_digits)


# Picks form with smallest |a|, positive if possible. If there is still
# a tie, pick smallest |b|.
def best_form(forms):
    best_a = -math.inf
    best_b = math.inf
    for (a, b, c) in forms:
        smaller_abs_a = (abs(a) < abs(best_a))
        same_abs_a_but_pos = ((abs(a) == abs(best_a)) and (a > best_a))
        same_a_smaller_abs_b = ((a == best_a) and (abs(b) < abs(best_b)))
        if smaller_abs_a or same_abs_a_but_pos or same_a_smaller_abs_b:
            best_a = a
            best_b = b
            best_c = c

    return (best_a, best_b, best_c)
 

# Find set of representatives of classes of quadratic forms of
# discriminant d, for d>0 not a square. For d fundamental,
# should give narrow class number of Q(sqrt(d)).
# See http://www.numbertheory.org/php/reduce.html
def reduced_forms(D):
    if (D % 4 != 0) and (D % 4 != 1):
        raise ValueError('D must be 0 or 1 mod 4', D)
    if gmpy2.is_square(D):
        raise ValueError('D cannot be a square', D)

    reduced_forms = []
    R = sqrt(D)
    for b in range(1, ceil(R)): # 0 < b < R
        for a in range(ceil((R-b)/2), ceil((R+b)/2)): # 0 < R - b < 2|a| < R + b
            if (b*b - D) % (4*a) == 0:
                c = (b*b - D) // (4*a)
                reduced_forms.append(((a, b, c), True))
                reduced_forms.append(((-a, b, -c), True))

    representatives = []
    for (form, is_new) in reduced_forms:
        if is_new:
            (equiv_forms, digits) = equiv_reduced_forms(form)
            best = best_form(equiv_forms) # pick best for aesthetics
            representatives.append((best, digits))
            for i in range(len(reduced_forms)):
                (f, f_is_new) = reduced_forms[i]
                if f_is_new and (f in equiv_forms):
                    reduced_forms[i] = (f, False)

    representatives.sort() # also for aesthetics
    return representatives


def print_forms(D):
    for (f, d) in reduced_forms(D):
        print(f)


def narrow_class_number(D):
    return len(reduced_forms(D))


def class_number(D):
    Pell_has_sol = False
    reps = reduced_forms(D)
    for (f, digits) in reps:
        if len(digits) % 2 == 1:
            Pell_has_sol = True

    if Pell_has_sol:
        return len(reps)
    else:
        return len(reps) // 2


def find_large_class_number(a):
    for t in range (a, a+10):
        d = 4*t*t + 1
        if is_square_free(d):
            print(d)
            print(narrow_class_number(d))
            print()

def class_number_search(a, b, n, narrow = False):
    res = []
    if a < 2:
        a=2
    for d in range(a, b):
        if d % 1000 == 0:
            print(d)
        if not is_square_free(d):
            continue
        if d % 4 == 2 or d % 4 == 3:
            if narrow:
                if narrow_class_number(4*d) == n:
                    res.append(d)
            else:
                if class_number(4*d) == n:
                    res.append(d)
        if d % 4 == 1:
            if narrow:
                if narrow_class_number(d) == n:
                    res.append(d)
            else:
                if class_number(d) == n:
                    res.append(d)

    return res


def make_arc(center, p, q, color = 'k', linewidth = 0.5):
    if q[0] > p[0]:
        p, q = q, p

    (px, py) = p
    (qx, qy) = q
        
    rp = sqrt((center - px)**2 + (py)**2)
    rq = sqrt((center - qx)**2 + (qy)**2)
    r = (rp + rq)/2
    if abs(rp - rq) > 10**(-10):
        raise ValueError('points do not form an arc')
    
    theta1 = np.rad2deg(np.arccos(min((px - center) / r, 1)))
    theta2 = np.rad2deg(np.arccos(max((qx - center) / r, -1)))

    return patches.Arc((center, 0), 2*r, 2*r, angle = 0, theta1 = theta1, theta2 = theta2,
                             linewidth = linewidth, color = color)


def intersects_left(form):
    (a, b, c) = form
    D = b*b - 4*a*c

    # center = -b / (2*a)
    # r_squared = D / (4*a*a)
    # check if (3/4 + (center + 1/2)**2 <= r_squared)
    # but use only integer operations
    
    return (3*a*a + (-b + a)*(-b + a) <= D)


def intersects_right(form):
    (a, b, c) = form
    D = b*b - 4*a*c

    # center = -b / (2*a)
    # r_squared = D / (4*a*a)
    # check if (3/4 + (center - 1/2)**2 <= r_squared)
    # but use only integer operations
    
    return (3*a*a + (b + a)*(b + a) <= D)


def intersects_bottom(form):
    (a, b, c) = form
    D = b*b - 4*a*c

    # center = -b / (2*a)
    # r_squared = D / (4*a*a)
    # check if one of (3/4 + (center + 1/2)**2 <= r_squared) or (3/4 + (center - 1/2)**2 <= r_squared) fails 
    # (with XOR as at least one is true) but use only integer operations
    
    return (3*a*a + (-b + a)*(-b + a) <= D) ^ (3*a*a + (b + a)*(b + a) <= D)


def has_intersection(form):
    return intersects_left(form) or intersects_right(form) or intersects_bottom(form)


def domain_intersection_pts(form):
    (a, b, c) = form
    D = b**2 - 4*a*c
    center = -b / (2*a)
    r_squared = D / (4*a*a)
    coords = []
    
    if intersects_left(form):
        px = -1/2
        py = sqrt(r_squared - (center + 1/2)**2)
        p = (px, py)
        coords.append(left)
    if intersects_bottom(form):
        rx = (1 + center**2 - r_squared) / (2*center)
        ry = sqrt(1 - rx**2)
        r = (rx, ry)
        coords.append(bottom)
        if intersects_left(form):
            q = r
        else:
            p = r
    if intersects_right(form):
        qx = 1/2
        qy = sqrt(r_squared - (center - 1/2)**2)
        q = (qx, qy)
        coords.append(right)

    return (p, q, coords)


def transform(form, last_int):
    (a, b, c) = form
    if last_int == left:
        res = (a, b - 2*a, a - b + c)
    elif last_int == bottom:
        res = (c, -b, a)
    else:
        res = (a, b + 2*a, a + b + c)

    return res


def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


def draw_geodesics(D):
    forms = reduced_forms(D)
    narrow_cn = len(forms)
    Pell_has_sol = False
    for (f, digits) in forms:
        if len(digits) % 2 == 1:
            Pell_has_sol = True
    if Pell_has_sol:
        cn = len(forms)
    else:
        cn = len(forms) // 2

    fig, ax = plt.subplots(figsize=(2,4))
##    plt.axis('off')
    ax.set_aspect("equal") # try "auto"
    ax.set_xlim((-2, 2))
    ax.set_ylim((0, 10))
##    ax.set(title=r'$D =$ %d: $h^+(D) =$ %d and $h(D) = $ %d' % (D, narrow_cn, cn))
    cmap = get_cmap(narrow_cn)
    
    # draw fundamental domain
    bottom_arc = make_arc(0, (-1/2, sqrt(3)/2), (1/2, sqrt(3)/2), 'k', 0.5)
    left_line = lines.Line2D([-1/2, -1/2], [sqrt(3)/2, 100], linewidth = 0.5, color = 'k')
    right_line = lines.Line2D([1/2, 1/2], [sqrt(3)/2, 100],  linewidth = 0.5, color = 'k')
    ax.add_patch(bottom_arc)
    ax.add_line(left_line)
    ax.add_line(right_line)

    # draw geodesics
    i = 0
    for (f, digits) in forms:
        color = cmap(i)
        i += 1
        
        form = f
        if not has_intersection(form):
                raise ValueError('Should not happen for our representatives, this is bad', form)
        if intersects_bottom(form):
            last_int = bottom
        else:
            if form[0] > 0:
                last_int = right
            else:
                last_int = left
        
        while True:
            if not has_intersection(form):
                raise ValueError('Should not happen for our representatives, this is bad', form)
            (a, b, c) = form
            center = -b / (2*a)
            (p, q, coords) = domain_intersection_pts(form)
            arc = make_arc(center, p, q, color, 0.3)
            ax.add_patch(arc)

##            root1 = (-b + sqrt(D)) / (2*a)
##            root2 = (-b - sqrt(D)) / (2*a)
##            arc = make_arc(center, (root1, 0), (root2, 0), 'b', 0.5)
##            ax.add_patch(arc)

            coords.remove(last_int)
            last_int = coords[0]
            form = transform(form, last_int)
            
            if last_int == bottom:
                last_int = bottom
            elif last_int == left:
                last_int = right
            elif last_int == right:
                last_int = left

##            plt.show(block=False)
##            plt.pause(0.1)

            if form == f:
                break

    # redraw fundamental domain
    bottom_arc = make_arc(0, (-1/2, sqrt(3)/2), (1/2, sqrt(3)/2), 'k', 0.5)
    left_line = lines.Line2D([-1/2, -1/2], [sqrt(3)/2, 100], linewidth = 0.5, color = 'k')
    right_line = lines.Line2D([1/2, 1/2], [sqrt(3)/2, 100],  linewidth = 0.5, color = 'k')
    ax.add_patch(bottom_arc)
    ax.add_line(left_line)
    ax.add_line(right_line)
    
    print(narrow_cn)

    fig.tight_layout()
    file_name = 'geodesics_' + str(D)
##    plt.savefig(file_name, format = 'pgf')
    plt.show()
