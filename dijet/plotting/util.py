# coding: utf-8

def eta_bin(eta_lo, eta_hi):
    if abs(eta_lo) < 1e-3:
        return fr"$|\eta|$ < {eta_hi:g}"
    else:
        return fr"{eta_lo:g} < $|\eta|$ < {eta_hi:g}"


def pt_bin(pt_lo, pt_hi):
    spt = r"$p_\mathrm{T}$"
    return fr"{int(pt_lo):d} < {spt} < {int(pt_hi):d} GeV"


def alpha_bin(alpha):
    salpha = r"$\alpha$ < "
    return fr"{salpha}{alpha}"


def dot_to_p(var):
    return f"{var}".replace(".", "p")


def add_text(ax, x, y, text, offset=0, ha="left", va="top"):
    ax.text(
        x, y - offset, text,
        transform=ax.transAxes,
        horizontalalignment=ha,
        verticalalignment=va,
    )
