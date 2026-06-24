

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

################################### plotting_rotation.##############################################

def function_to_latex(function_plot):
    """
    Convert a function name into a clean LaTeX label.
    """

    if function_plot == "constant_0":
        return r"\mathrm{constant}_0"

    if function_plot == "constant_1":
        return r"\mathrm{constant}_1"

    if function_plot == "balanced":
        return r"\mathrm{balanced}"

    return rf"\mathrm{{{function_plot}}}"


def error_label(error_pos, axis):
    """
    Convert an error-position name into a readable plot label.
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
    """
    Add the ideal-circuit probability P0_ideal to each matching row.

    This lets us compare noisy/error results against the no_error circuit.
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

def rotation_error_effect(df, n_plot, axis_plot, function_plot, shots):
    """
    Plot P(0...0) against rotation angle for different error positions.
    """

    # Keep only the data needed for this plot.
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

    # Plot one curve for each error position.
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
        f"{axis}-axis rotation error in the Deutsch--Jozsa algorithm\n"
        f"function = {function_plot}, n = {n_plot}, shots = {shots}",
        fontsize=13
    )

    plt.xlim(0, 180)
    plt.ylim(-0.1, 1.05)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()

def rotation_error_sensitivity(
    df,
    n_plot,
    axis_plot,
    target_plot,
    function_plot,
    shots
):
    """
    Plot sensitivity to rotation errors using success probability.

    Sensitivity = |success_error - success_ideal|
    """

    axis = str(axis_plot).strip().upper()

    ideal = (
        df[df["error_position"] == "no_error"]
        [["n", "theta_deg", "axis", "target_qubit", "function", "success"]]
        .rename(columns={"success": "success_ideal"})
    )

    df2 = df.merge(
        ideal,
        on=["n", "theta_deg", "axis", "target_qubit", "function"],
        how="left"
    )

    plot_df = df2[
        (df2["n"] == n_plot)
        & (df2["axis"] == axis)
        & (df2["target_qubit"] == target_plot)
        & (df2["function"] == function_plot)
        & (df2["error_position"] != "no_error")
    ].copy()

    if plot_df.empty:
        print("No matching data found.")
        return

    plot_df["delta_success"] = np.abs(
        plot_df["success"] - plot_df["success_ideal"]
    )

    plt.figure(figsize=(8, 5), dpi=150)

    for error_pos in [
        "E1_before_H",
        "E2_after_first_H",
        "E3_after_oracle",
        "E4_after_final_H",
    ]:
        subset = plot_df[
            plot_df["error_position"] == error_pos
        ].sort_values("theta_deg")

        if subset.empty:
            continue

        plt.plot(
            subset["theta_deg"],
            subset["delta_success"],
            marker="o",
            linewidth=2,
            markersize=5,
            label=error_label(error_pos, axis),
        )

    function_label = function_to_latex(function_plot)

    plt.xlabel(r"Rotation angle $\theta$ (degrees)", fontsize=12)
    plt.ylabel(
        r"Sensitivity $|S_{\mathrm{error}} - S_{\mathrm{ideal}}|$",
        fontsize=12,
    )

    plt.title(
        rf"Sensitivity to $R_{{{axis}}}(\theta)$ errors"
        "\n"
        rf"$f(x)={function_label},\; n={n_plot},\; "
        rf"q={target_plot},\; \mathrm{{shots}}={shots}$",
        fontsize=13,
    )

    plt.xlim(0, 180)
    plt.ylim(-0.05, 1.05)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()

def rotation_axis_comparison(df, n_plot, function_plot, error_position, shots):
    """
    Compare X, Y, and Z rotation errors at one error position.
    """

    plot_df = df[
        (df["n"] == n_plot) &
        (df["function"] == function_plot) &
        (df["error_position"] == error_position)
    ]

    plt.figure(figsize=(8, 5), dpi=150)

    # Plot one curve for each rotation axis.
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

    # Fixed broken f-string here.
    plt.title(
        rf"Comparison of $R_X$, $R_Y$, and $R_Z$ errors at {error_position}"
        "\n"
        rf"$f(x)={function_label},\; n={n_plot},\; \mathrm{{shots}}={shots}$",
        fontsize=13
    )

    plt.xlim(0, 180)
    plt.ylim(-0.1, 1.05)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show() 


def rotation_error_scalability(df, theta_choice, axis_choice, error_choice, function_choice, shots):
    """
    Plot how success probability changes as n increases for one function.
    """

    axis = str(axis_choice).strip().upper()

    labels = {
        "no_error": "Ideal circuit",
        "E1_before_H": rf"$R_{{{axis}}}(\theta)$ at $E_1$",
        "E2_after_first_H": rf"$R_{{{axis}}}(\theta)$ at $E_2$",
        "E3_after_oracle": rf"$R_{{{axis}}}(\theta)$ at $E_3$",
        "E4_after_final_H": rf"$R_{{{axis}}}(\theta)$ at $E_4$",
    }

    scale_df = df[
        (df["theta_deg"] == theta_choice)
        & (df["axis"] == axis)
        & (df["error_position"] == error_choice)
        & (df["function"] == function_choice)
    ]

    grouped = (
        scale_df
        .groupby("n", as_index=True)["success"]
        .mean()
        .sort_index()
    )

    plt.figure(figsize=(8, 5), dpi=150)

    plt.plot(
        grouped.index,
        grouped.values,
        marker="o",
        linewidth=2,
        markersize=5,
    )

    plt.xlabel("Number of input qubits n", fontsize=12)
    plt.ylabel("Success probability", fontsize=12)

    plt.title(
    f"Scalability for {function_choice} under {labels.get(error_choice, error_choice)} noise\n"
    f"axis={axis}, θ={theta_choice}°, shots={shots}",
    fontsize=13,)

    plt.ylim(0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    filename = f"scalability_{function_choice}_{error_choice}_{axis}_{theta_choice}deg.pdf"

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.show()


def rotation_error_success_probability(
    df,
    n_plot,
    axis_plot,
    target_plot,
    function_plot,
    shots
):
    """
    Plot success probability versus rotation angle for one Deutsch–Jozsa function.

    Parameters
    ----------
    df : pandas.DataFrame
        Results dataframe containing:
        n, theta_deg, axis, target_qubit,
        function, error_position, success, shots.

    n_plot : int
        Number of input qubits.

    axis_plot : str
        Rotation axis ("X", "Y", or "Z").

    target_plot : int
        Target qubit where the rotation error is applied.

    function_plot : str
        One of:
        "constant_0", "constant_1", "balanced".

    shots : int
        Number of measurement shots.
    """

    axis = axis_plot.upper()

    plot_df = df[
        (df["n"] == n_plot)
        & (df["axis"] == axis)
        & (df["target_qubit"] == target_plot)
        & (df["function"] == function_plot)
    ]

    if plot_df.empty:
        print("No matching data found.")
        return

    labels = {
        "no_error": "Ideal circuit",
        "E1_before_H": rf"$R_{{{axis}}}(\theta)$ at $E_1$",
        "E2_after_first_H": rf"$R_{{{axis}}}(\theta)$ at $E_2$",
        "E3_after_oracle": rf"$R_{{{axis}}}(\theta)$ at $E_3$",
        "E4_after_final_H": rf"$R_{{{axis}}}(\theta)$ at $E_4$",
    }

    plt.figure(figsize=(10, 6), dpi=150)

    for error_pos, label in labels.items():

        subset = (
            plot_df[plot_df["error_position"] == error_pos]
            .sort_values("theta_deg")
        )

        if subset.empty:
            continue

        plt.plot(
            subset["theta_deg"],
            subset["success"],
            marker="o",
            markersize=5,
            linewidth=2,
            label=label,
        )

    plt.xlabel(r"Rotation angle $\theta$ (degrees)", fontsize=12)
    plt.ylabel("P(0...0)", fontsize=12)

    plt.title(
        f"Probability under {axis}-Axis Rotation Errors\n"
        f"{function_plot.replace('_', ' ').title()}, "
        f"n={n_plot}, target qubit={target_plot}, shots={shots}",
        fontsize=13,
    )

    plt.xlim(0, 180)
    plt.ylim(0, 1.05)

    plt.grid(True, alpha=0.3)
    plt.legend(title="Error position")
    plt.tight_layout()

    filename = (
        f"rotation_success_"
        f"{function_plot}_"
        f"n{n_plot}_"
        f"{axis}_"
        f"q{target_plot}.pdf"
    )

    plt.savefig(
        filename,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


############################plotting_depolarizing###############################

def depolarizing_noise_effect(df, n_plot, function_plot, shots=None):
    """
    Plot success probability against depolarizing probability p.
    """

    plot_df = df[
        (df["n"] == n_plot) &
        (df["function"] == function_plot)
    ]

    plt.figure(figsize=(8, 5), dpi=150)

    labels = {
        "no_error": "Ideal circuit",
        "D(P)1_before_H": r"$D_p$ at $E_1$ before H",
        "D(P)2_after_first_H": r"$D_p$ at $E_2$ after first H",
        "D(P)3_after_oracle": r"$D_p$ at $E_3$ after oracle",
        "D(P)4_after_final_H": r"$D_p$ at $E_4$ after final H",
    }

    for error_pos, label in labels.items():
        subset = plot_df[plot_df["error_position"] == error_pos].sort_values("p")

        if subset.empty:
            continue

        plt.plot(
            subset["p"],
            subset["success"],
            marker="o",
            linewidth=2,
            markersize=5,
            label=label
        )

    function_label = function_to_latex(function_plot)
    title_shots = f", shots={shots}" if shots is not None else ""

    plt.xlabel(r"Depolarizing probability $p$", fontsize=12)
    plt.ylabel("Success probability", fontsize=12)

    plt.title(
        r"Effect of depolarizing noise on the Deutsch--Jozsa algorithm"
        "\n"
        rf"$f(x)={function_label},\; n={n_plot}{title_shots}$",
        fontsize=13
    )

    plt.xlim(0, 1.0)
    plt.ylim(-0.05, 1.05)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()

def depolarizing_noise_sensitivity(df, n_plot, function_plot, shots):
    """
    Plot sensitivity to depolarizing noise.

    Sensitivity means:

        |success_noise - success_ideal|
    """

    ideal = (
        df[df["error_position"] == "no_error"]
        [["n", "p", "function", "success"]]
        .rename(columns={"success": "success_ideal"})
    )

    df2 = df.merge(
        ideal,
        on=["n", "p", "function"],
        how="left"
    )

    plot_df = df2[
        (df2["n"] == n_plot) &
        (df2["function"] == function_plot) &
        (df2["error_position"] != "no_error")
    ].copy()

    plot_df["delta_success"] = np.abs(
        plot_df["success"] - plot_df["success_ideal"]
    )

    plt.figure(figsize=(8, 5), dpi=150)

    labels = {
        "D(P)1_before_H": r"$D_p$ at $E_1$ before H",
        "D(P)2_after_first_H": r"$D_p$ at $E_2$ after first H",
        "D(P)3_after_oracle": r"$D_p$ at $E_3$ after oracle",
        "D(P)4_after_final_H": r"$D_p$ at $E_4$ after final H",
    }

    for error_pos, label in labels.items():
        subset = plot_df[plot_df["error_position"] == error_pos].sort_values("p")

        if subset.empty:
            continue

        plt.plot(
            subset["p"],
            subset["delta_success"],
            marker="o",
            linewidth=2,
            markersize=5,
            label=label
        )

    function_label = function_to_latex(function_plot)

    plt.xlabel(r"Depolarizing probability $p$", fontsize=12)
    plt.ylabel(r"Sensitivity $|S_{\mathrm{noise}} - S_{\mathrm{ideal}}|$", fontsize=12)

    plt.title(
        r"Sensitivity to depolarizing noise in the Deutsch--Jozsa algorithm"
        "\n"
        rf"$f(x)={function_label},\; n={n_plot},\; \mathrm{{shots}}={shots}$",
        fontsize=13
    )

    plt.xlim(0, 1.0)
    plt.ylim(-0.05, 1.05)
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()