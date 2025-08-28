import streamlit as st
import sympy as smp
import numpy as np
import matplotlib.pyplot as plt

t = smp.symbols('t', real=True)

# ---- from your script: helper to build the homogeneous solution Xi ----
def solver(a1, b1, c1):
    if a1 == 0:
        # 1st-order case
        return smp.exp(-c1*t/b1)
    # 2nd-order case
    disc = b1**2 - 4*a1*c1
    x1 = (-b1 + smp.sqrt(disc)) / (2*a1)
    x2 = (-b1 - smp.sqrt(disc)) / (2*a1)
    if smp.simplify(x1 - x2) == 0:
        m, n = 0, 1
        return (m + n*t) * smp.exp(x1*t)
    else:
        m = -1/(x2 - x1)
        n =  1/(x2 - x1)
        return m*smp.exp(x1*t) + n*smp.exp(x2*t)

def unit_impulse_response(a1, b1, c1, a2, b2, c2, symbolic_unit=True):
    # h(t) = a2*δ(t) + (a2*x'' + b2*x' + c2*x)*u(t)  (2nd-order)
    #       b2*δ(t) + (a2*x'' + b2*x' + c2*x)*u(t)    (1st-order, a1=0)
    Xi = solver(a1, b1, c1)
    u = smp.Function('u')(t) if symbolic_unit else 1
    δ = smp.Symbol('δ(t)') if symbolic_unit else 0
    if a1 == 0:
        return b2*δ + (a2*smp.diff(Xi, t, 2) + b2*smp.diff(Xi, t) + c2*Xi)*u
    else:
        return a2*δ + (a2*smp.diff(Xi, t, 2) + b2*smp.diff(Xi, t) + c2*Xi)*u

def zero_state_response(a1, b1, c1, a2, b2, c2):
    Xi = solver(a1, b1, c1)
    integrand = a2*smp.diff(Xi, t, 2) + b2*smp.diff(Xi, t) + c2*Xi
    if a1 == 0:
        return b2 + smp.integrate(integrand, (t, 0, t))
    else:
        return a2 + smp.integrate(integrand, (t, 0, t))

def _plot_expr(expr, T=50, title="", ylabel=""):
    # Plot a sympy expression by sampling to numpy and using matplotlib
    f = smp.lambdify(t, smp.simplify(expr), "numpy")
    ts = np.linspace(0, T, 1000)
    ys = f(ts)
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot(ts, ys)
    ax.set_title(title)
    ax.set_xlabel("Time t")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    st.pyplot(fig)
    plt.close(fig)

st.set_page_config(page_title="LTIC Responses", page_icon="⚙️", layout="centered")
st.title("LTIC Response Visualizer")
st.caption("Compute Zero-Input, Impulse, Step, and Total Output for first/second-order systems.")

with st.form("coeffs"):
    st.subheader("Enter coefficients")
    col1, col2 = st.columns(2)
    with col1:
        a1 = st.number_input("a (coeff of D² in Q(D))", value=1.0)
        b1 = st.number_input("b (coeff of D in Q(D))",  value=2.0)
        c1 = st.number_input("c (constant in Q(D))",     value=1.0)
    with col2:
        a2 = st.number_input("d (coeff of D² in P(D))",  value=0.0)
        b2 = st.number_input("e (coeff of D in P(D))",   value=1.0)
        c2 = st.number_input("f (constant in P(D))",     value=0.0)
    T = st.slider("Plot duration T (seconds)", 10, 200, 50, 5)
    submitted = st.form_submit_button("Compute")

if submitted:
    # basic sanity checks (match your script’s constraints)
    if (a1 == 0 and a2 != 0) or (a1 == b1 == 0) or (a2 == b2 == 0):
        st.error("Please check orders: require P(D) order ≤ Q(D) order; min/max orders of Q(D), P(D) are 1 and 2.")
    else:
        Xi = solver(a1, b1, c1)
        h_sym = unit_impulse_response(a1, b1, c1, a2, b2, c2, symbolic_unit=True)
        h_num = unit_impulse_response(a1, b1, c1, a2, b2, c2, symbolic_unit=False)
        Xs = zero_state_response(a1, b1, c1, a2, b2, c2)
        Xo = smp.simplify(Xi + Xs)

        st.markdown("### Symbolic results")
        st.latex(r"\text{Zero-input } \Xi(t) = " + smp.latex(smp.simplify(Xi)))
        st.latex(r"\text{Impulse } h(t) = " + smp.latex(h_sym))
        st.latex(r"\text{Step } X_{zs}(t) = " + smp.latex(smp.simplify(Xs)))
        st.latex(r"\text{Output } X(t) = " + smp.latex(Xo) + r"\;\cdot u(t)")

        st.markdown("### Plots")
        _plot_expr(Xo, T=T, title="Output X(t)", ylabel="x(t)")
        _plot_expr(h_num, T=T, title="Impulse Response h(t)", ylabel="h(t)")
        _plot_expr(Xs, T=T, title="Step Response Xzs(t)", ylabel="Xzs(t)")
