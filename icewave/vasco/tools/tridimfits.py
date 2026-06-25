#%%
import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#%%

## ── Données (remplacez par les vôtres) ────────────────────────────────────────
#rng = np.random.default_rng(42)
#x = rng.uniform(-3, 3, 80)
#y = rng.uniform(-3, 3, 80)
#z = 1.5*x**2 - 0.8*y**2 + 0.4*x*y + 2*x - y + 3 + rng.normal(0, 0.5, 80)

def fit_3dparabola(x,y,z, plot=False):

    XY = np.column_stack([x, y])   # shape (n, 2)

    # ── Fit en 2 lignes ───────────────────────────────────────────────────────────
    modele = make_pipeline(PolynomialFeatures(degree=2), LinearRegression())
    modele.fit(XY, z)

    # ── Résultats ─────────────────────────────────────────────────────────────────
    z_pred = modele.predict(XY)
    print(f"R²   : {r2_score(z, z_pred):.4f}")
    print(f"RMSE : {np.sqrt(np.mean((z - z_pred)**2)):.4f}")
    print(f"Coefficients : {modele[-1].coef_}")
    print(f"Intercept    : {modele[-1].intercept_:.4f}")

    # ── Visualisation ─────────────────────────────────────────────────────────────
    xg, yg = np.meshgrid(np.linspace(np.min(x), np.max(x), 50), np.linspace(np.min(y), np.max(y), 50))
    Zg = modele.predict(np.column_stack([xg.ravel(), yg.ravel()])).reshape(xg.shape)
    
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(xg, yg, Zg, alpha=0.4, cmap='viridis')
        ax.scatter(x, y, z, c='red', s=20)
        
        plt.show()
    
    dict_results = {}
    dict_results['coefficients'] = modele[-1].coef_
    dict_results['intercept'] = modele[-1].intercept_
    dict_results['r2_score'] = r2_score(z, z_pred)
    dict_results['RMS'] = np.sqrt(np.mean((z - z_pred)**2))
    dict_results['order_coef'] = ['1','x','y','x²','xy','y²']
    dict_results['info'] = 'intercept is not coef order 1 here but in intercept key'

    return dict_results


def measure_bidimensional_curvature_around_central_point(matrix2d, xind_centre, yind_centre, window_size=(11,11)):

    dx = int((window_size[0]-1)/2)
    dy = int((window_size[1]-1)/2)


    xind_inf = xind_centre-dx
    xind_sup = xind_centre+dx
    yind_inf = yind_centre-dy
    yind_sup = yind_centre+dy


    xvals = np.arange(uz.shape[1])
    yvals = np.arange(uz.shape[0])

    X,Y = np.meshgrid(xvals,yvals)
    X = X[yind_inf:yind_sup, xind_inf:xind_sup]
    Y = Y[yind_inf:yind_sup, xind_inf:xind_sup]

    x = X.flatten()
    y = Y.flatten()

    Z = matrix2d[yind_inf:yind_sup,xind_inf:xind_sup]
    z = Z.flatten()

    dict_results = fit_3dparabola(x,y,z, plot=True)
    return dict_results

import numpy as np
from sympy import symbols, cos, sin, expand, collect, Symbol

def convertir_base_tournee(coeffs, alpha):
    """
    coeffs : [a, b, c, d, e, f] du modèle
             f(x,y) = a*x² + b*y² + c*x*y + d*x + e*y + f
    alpha  : angle de rotation (radians)
    
    Retourne les nouveaux coefficients [A, B, C, D, E, F]
    dans la base tournée (u, v).
    """
    a, b, c, d, e, f = coeffs
    u, v = symbols('u v')
    ca, sa = np.cos(alpha), np.sin(alpha)

    # Substitution : x et y en fonction de u, v
    x_expr = u * ca - v * sa
    y_expr = u * sa + v * ca

    # Reconstruction du polynôme dans la nouvelle base
    poly = (a * x_expr**2
          + b * y_expr**2
          + c * x_expr * y_expr
          + d * x_expr
          + e * y_expr
          + f)

    poly = expand(poly)

    # Extraction des coefficients
    A = float(poly.coeff(u, 2))          # u²
    B = float(poly.coeff(v, 2))          # v²
    C = float(poly.coeff(u).coeff(v))    # u·v
    D = float(poly.coeff(u, 1).subs(v, 0))  # u
    E = float(poly.coeff(v, 1).subs(u, 0))  # v
    F = float(poly.subs({u: 0, v: 0}))   # constante

    return np.array([A, B, C, D, E, F])



def change_order_from_polyfeatures_to_polycoefs(coef_polyfeatures=np.zeros(6), intercept=0):
    a = coef_polyfeatures[3]
    b = coef_polyfeatures[5]
    c = coef_polyfeatures[4]
    d = coef_polyfeatures[1]
    e = coef_polyfeatures[2]
    if intercept==None:
        f = coef_polyfeatures[0]
    else:
        f = intercept
    coefficients = np.array([a,b,c,d,e,f])
    return coefficients
# %%
