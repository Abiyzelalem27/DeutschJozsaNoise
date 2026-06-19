


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 


def plot_error_positions(df, n_plot, axis_plot, function_plot):

    plot_df = df[
        (df["n"] == n_plot) &
        (df["axis"] == axis_plot) &
        (df["function"] == function_plot)
    ]

    plt.figure(figsize=(8, 5), dpi=150)

    axis = axis_plot.upper()

    labels = {
        "no_error": "Ideal circuit",
        "E1_before_H": rf"$R_{{{axis}}}\!\left(\frac{{\pi}}{{2}}\right)$ at $E_1$",
        "E2_after_first_H": rf"$R_{{{axis}}}\!\left(\frac{{\pi}}{{2}}\right)$ at $E_2$",
        "E3_after_oracle": rf"$R_{{{axis}}}\!\left(\frac{{\pi}}{{2}}\right)$ at $E_3$",
        "E4_after_final_H": rf"$R_{{{axis}}}\!\left(\frac{{\pi}}{{2}}\right)$ at $E_4$",
    }

    for error_pos, label in labels.items():

        subset = plot_df[plot_df["error_position"] == error_pos]

        if subset.empty:
            continue

        plt.plot(
            subset["theta_deg"],
            subset["P0"],
            marker="o",
            linewidth=2,
            markersize=5,
            label=label
        )

    
    plt.xlabel(r"Rotation angle $\theta$ (degrees)", fontsize=12)
    plt.ylabel(r"Probability $P(0\ldots0)$", fontsize=12)
    plt.title(
    f"{axis_plot}-axis rotation error in the Deutsch--Jozsa algorithm\n"
    f"parity = {function_plot}, n = {n_plot}, shots = {shots}",
    fontsize=13,)

    plt.ylim(-0.1, 1.05)
    plt.xlim(0, 180)

    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()



def function_to_latex(function_plot):
    """
    Convert function names to LaTeX labels for figure titles.

    Parameters
    ----------
    function_plot : str
        Function type.

    Returns
    -------
    str
        LaTeX-formatted function label.
    """

    if function_plot == "constant_0":
        return r"\mathrm{constant}_0"

    elif function_plot == "constant_1":
        return r"\mathrm{constant}_1"

    elif function_plot == "balanced":
        return r"\mathrm{balanced}"

    else:
        return rf"\mathrm{{{function_plot}}}"


def error_label(error_pos, axis):
    """
    Convert error position names to publication-quality legend labels.

    Parameters
    ----------
    error_pos : str
        Error location identifier.
    axis : str
        Rotation axis ("X", "Y", or "Z").

    Returns
    -------
    str
        LaTeX-formatted legend label.
    """

    if error_pos == "no_error":
        return "Ideal circuit"

    position_map = {
        "E1_before_H": "E_1",
        "E2_after_first_H": "E_2",
        "E3_after_oracle": "E_3",
        "E4_after_final_H": "E_4",
    }

    return (
        rf"$R_{{{axis}}}\!\left(\frac{{\pi}}{{2}}\right)$ "
        rf"at ${position_map[error_pos]}$"
    )


def add_ideal_reference(df):
    ideal = df[df["error_position"] == "no_error"][
        ["n", "theta_deg", "axis", "function", "P0"]
    ].rename(columns={"P0": "P0_ideal"})

    return df.merge(
        ideal,
        on=["n", "theta_deg", "axis", "function"],
        how="left"
    )


def plot_error_positions(df, n_plot, axis_plot, function_plot, shots=1024):

    plot_df = df[
        (df["n"] == n_plot) &
        (df["axis"] == axis_plot) &
        (df["function"] == function_plot)
    ]

    plt.figure(figsize=(8, 5), dpi=150)

    axis = axis_plot.upper()

    labels = {
        "no_error": "Ideal circuit",
        "E1_before_H": rf"$R_{{{axis}}}\!\left(\frac{{\pi}}{{2}}\right)$ at $E_1$",
        "E2_after_first_H": rf"$R_{{{axis}}}\!\left(\frac{{\pi}}{{2}}\right)$ at $E_2$",
        "E3_after_oracle": rf"$R_{{{axis}}}\!\left(\frac{{\pi}}{{2}}\right)$ at $E_3$",
        "E4_after_final_H": rf"$R_{{{axis}}}\!\left(\frac{{\pi}}{{2}}\right)$ at $E_4$",
    }

    for error_pos, label in labels.items():
        subset = plot_df[plot_df["error_position"] == error_pos]

        if subset.empty:
            continue

        plt.plot(
            subset["theta_deg"],
            subset["P0"],
            marker="o",
            linewidth=2,
            markersize=5,
            label=label
        )

    function_label = function_to_latex(function_plot)

    plt.xlabel(r"Rotation angle $\theta$ (degrees)", fontsize=12)
    plt.ylabel(r"Probability $P(0\ldots0)$", fontsize=12)

    plt.title(
        rf"{axis_plot}-axis rotation error in the Deutsch--Jozsa algorithm"
        "\n"
        rf"$f(x)={function_label},\; n={n_plot},\; \mathrm{{shots}}={shots}$",
        fontsize=13,
    )

    plt.ylim(-0.1, 1.05)
    plt.xlim(0, 180)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()


def plot_sensitivity(df, n_plot, axis_plot, function_plot, shots=1024):
    """
    Plot sensitivity of Deutsch-Jozsa success probability to rotation errors.

    Sensitivity is defined as:

        |P_error - P_ideal|

    where P_ideal is the probability from the no_error circuit.

    Parameters
    ----------
    df : pandas.DataFrame
        Simulation results dataframe.
    n_plot : int
        Number of input qubits.
    axis_plot : str
        Rotation axis: "X", "Y", or "Z".
    function_plot : str
        Function type: "constant_0", "constant_1", or "balanced".
    shots : int
        Number of measurement shots.
    """

    df2 = add_ideal_reference(df)

    plot_df = df2[
        (df2["n"] == n_plot) &
        (df2["axis"] == axis_plot) &
        (df2["function"] == function_plot) &
        (df2["error_position"] != "no_error")
    ].copy()

    plot_df["delta_P"] = np.abs(plot_df["P0"] - plot_df["P0_ideal"])

    plt.figure(figsize=(8, 5), dpi=150)

    axis = axis_plot.upper()

    for error_pos in [
        "E1_before_H",
        "E2_after_first_H",
        "E3_after_oracle",
        "E4_after_final_H",
    ]:
        subset = plot_df[plot_df["error_position"] == error_pos]

        if subset.empty:
            continue

        plt.plot(
            subset["theta_deg"],
            subset["delta_P"],
            marker="o",
            linewidth=2,
            markersize=5,
            label=error_label(error_pos, axis)
        )

    function_label = function_to_latex(function_plot)

    plt.xlabel(r"Rotation angle $\theta$ (degrees)", fontsize=12)
    plt.ylabel(r"Sensitivity $|P_{\mathrm{error}} - P_{\mathrm{ideal}}|$", fontsize=12)

    plt.title(
        rf"Sensitivity to $R_{{{axis}}}(\theta)$ errors"
        "\n"
        rf"$f(x)={function_label},\; n={n_plot},\; \mathrm{{shots}}={shots}$",
        fontsize=13,
    )

    plt.xlim(0, 180)
    plt.ylim(-0.05, 1.05)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()


def plot_axis_comparison(df, n_plot, function_plot, error_position, shots=1024):

    plot_df = df[
        (df["n"] == n_plot) &
        (df["function"] == function_plot) &
        (df["error_position"] == error_position)
    ]

    plt.figure(figsize=(8, 5), dpi=150)

    for axis_name in ["X", "Y", "Z"]:
        subset = plot_df[plot_df["axis"] == axis_name]

        if subset.empty:
            continue

        plt.plot(
            subset["theta_deg"],
            subset["P0"],
            marker="o",
            linewidth=2,
            markersize=5,
            label=rf"$R_{{{axis_name}}}(\theta)$"
        )

    function_label = function_to_latex(function_plot)

    plt.xlabel(r"Rotation angle $\theta$ (degrees)", fontsize=12)
    plt.ylabel(r"Probability $P(0\ldots0)$", fontsize=12)

    plt.title(
        rf"Comparison of $R_X$, $R_Y$, and $R_Z$ errors at {error_position}"
        "\n"
        rf"$f(x)={function_label},\; n={n_plot},\; \mathrm{{shots}}={shots}$",
        fontsize=13,
    )

    plt.xlim(0, 180)
    plt.ylim(-0.1, 1.05)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()



def add_ideal_reference(df):
    """
    Add the ideal-circuit probability to every row.

    Creates a new column:
        P0_ideal
    """

    ideal = (
        df[df["error_position"] == "no_error"]
        [["n", "theta_deg", "axis", "function", "P0"]]
        .rename(columns={"P0": "P0_ideal"})
    )

    return df.merge(
        ideal,
        on=["n", "theta_deg", "axis", "function"],
        how="left"
    ) 

def plot_depolarizing_positions(df, n_plot, function_plot):

    plot_df = df[
        (df["n"] == n_plot) &
        (df["function"] == function_plot)
    ]

    plt.figure(figsize=(8, 5), dpi=150)

    labels = {
    "no_error": "Ideal circuit",
    "D(P)1_before_H": r"$D_p$ at $E_1$ (before H)",
    "D(P)2_after_first_H": r"$D_p$ at $E_2$ (after first H)",
    "D(P)3_after_oracle": r"$D_p$ at $E_3$ (after oracle)",
    "D(P)4_after_final_H": r"$D_p$ at $E_4$ (after final H)",
    }

    for error_pos, label in labels.items():
        subset = plot_df[plot_df["error_position"] == error_pos]

        plt.plot(
            subset["p"],
            subset["P0"],
            marker="o",
            linewidth=2,
            markersize=5,
            label=label
        )

    plt.xlabel(r"Depolarizing probability $p$", fontsize=12)
    plt.ylabel(r"Probability $P(0\ldots0)$", fontsize=12)
    plt.title(rf"Effect of depolarizing noise on the Deutsch--Jozsa algorithm"
    "\n"
    rf"$\mathrm{{function}}={function_plot},\ n={n_plot},\ "
    rf"\mathrm{{shots}}={4096}$",
    fontsize=13)
    
    plt.ylim(-0.02, 1.05)
    plt.xlim(0, 1.0)

    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show() 
