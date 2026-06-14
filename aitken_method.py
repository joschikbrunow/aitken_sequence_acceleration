from math import *
import matplotlib.pyplot as plt

# Class implementation of the Aitken-Delta-Square-Method
# Usable for explicit and recursive functions
class Sequence_Aitken:
    def __init__(self, A, lim, min_m, max_m, rec, x_0=0.0, autocreate=True):
        """
        Class for the usage of the explicit and recursive Aitken-Delta-Square-Method.

        Input:
            - A[sequence]: Sequence whose convergence is to be accelerated
            - lim[real]: Suspected limit of A
            - min_m[integer]: Index of the first term of A
            - max_m[integer]: Index of the last calculated term of A
            - rec[boolean]: True, if sequence is recursive
            - x_0[real/Optional]: Initial value in the case of a recursive function
            - autocreate[boolean/Optional]: True, if sequence is supposed to be calculated from index min_m to max_m at initilization
        """
        # Sequence formula/function
        self.A = A
        # Suspected limit
        self.lim = lim
        # Correct boundary indices if min is larger than max
        if (min_m <= max_m):
            self.min_m = min_m
            self.max_m = max_m
        else:
            self.min_m = max_m
            self.max_m = min_m
        # Space for calculated sequence
        self.sequence = []
        # Space for Aitken's sub-sequences (normal)
        self.aitken_sequence_normal = []
        # Space for Aitken's sub-sequences (recursive fixpoint iteration)
        self.aitken_sequence_fixedpoint = []
        # Bool if recursive
        self.rec = rec
        # Initial value in case of a recursive sequence
        self.x_0 = x_0

        # Create sequence at initilization if wanted
        if (autocreate):
            if (self.rec):
                self.set_sequence_recursive(self.x_0)
            else:
                self.set_sequence_normal()

    # Print important information of the class
    def __str__(self):
        """
        Prints important information of the class.
        
        Output:
            - output[String]: Information of the class
        """
        output = ""
        # Standard values
        output += f"Suspected limit {self.lim}\n"
        output += f"Minimal index: {self.min_m}\n"
        output += f"Maximal index: {self.max_m}\n"

        # Overview of the sequence
        if (len(self.sequence) != 0):
            s = "First values of the sequence: "
            for i in range(min(3, len(self.sequence))):
                s += f"{self.sequence[i]}, "
            
            output +=  f"{s[:-2]}\n"


            s = "Last values of the sequence: "
            for i in range(min(3, len(self.sequence))):
                    s += f"{self.sequence[len(self.sequence) - 3 + i]}, "

            output += f"{s[:-2]}\n"

        # Best computed limits
        try: # Limit already calculated (normal)
            output += f"Letzter Wert Aitken Normal: {self.aitken_sequence_normal[-1][-1]}\n"
        except: # Nothing to show
            pass
        try: # Limit already calculated (recursive)
            output += f"Letzter Wert Aitken Fixpunktbeschleunigt: {self.aitken_sequence_fixedpoint[-1]}\n"
        except: # Nothing to show
            pass
        return(output)

    # Creates sequence normal
    def set_sequence_normal(self):
        """
        Creates the explicit sequence by calculating formula for all indices.
        
        """
        for m in range(self.min_m, self.max_m+1):
            self.sequence.append(self.A(m))

    # Creates sequence recursive
    def set_sequence_recursive(self, x_0=None):
        """
        Creates the recursive sequence by calculating formula for the given range of indices.
        
        """
        # Use x_0 from class if it's not further specified
        if x_0 is None:
            x_0 = self.x_0

        # Initialize sequence with first value
        self.sequence = [x_0]

        # Calculate new term of the sequence with the term calculated before
        for _ in range(self.min_m, self.max_m):
            self.sequence.append(self.A(self.sequence[-1]))

    # Delta-operator
    def delta_A(self, A, m):
        """
        Delta-operator A_{m+1} - A_{m} for the given sequence A at the index m
        
        Input:
            - A[arr]: Sequence whose delta is calculated
            - m[integer]: Position where delta is computed
        Output:
            - Delta[real]: A_{m+1} - A_{m}
        """
        return(A[m+1] - A[m])
    
    # Delta-square-operator
    def delta_2_A(self, A, m):
        """
        Delta-square-operator A_{m+2} - 2 A_{m+1} + A_{m} for the given sequence A at the index m
        
        Input:
            - A[arr]: Sequence whose delta is calculated
            - m[integer]: Position where delta is computed
        Output:
            - Delta-square[real]: A_{m+2} - 2 A_{m+1} + A_{m}
        """
        return (self.delta_A(A, m+1) - self.delta_A(A, m))
    
    # Aitken-Delta-Square formula
    def aitken(self, A, m):
        """
        Aitken-method for sequence A at position m.

        A_m - (delta(A_m)^2)/(delta^2(A_m))
        
        Input:
            - A[arr]: Sequence whose delta is calculated
            - m[integer]: Position where delta is computed
        Output:
            - Aitken[real]: A_m - (delta(A_m)^2)/(delta^2(A_m))
        """
        d1 = self.delta_A(A, m)
        d2 = self.delta_2_A(A, m)
        try:
            return A[m] - (d1**2) / d2
        except ZeroDivisionError:
            return A[m]
    
    # Computes full explicit Aitken-sequence of A
    def aitken_normal(self, A):
        """
        Computes full explicit Aitken-sub-sequence of A.
        
        Input:
            - A[arr]: Sequence whose Aitken-sub-sequence is calculated
        Output:
            - Aitken-sub-sequence[arr]: Mit Aitken-Delta-Quadrat-Verfahren transformierte Folge
        """
        aitken_sequence = []
        # Using the Aitken scheme on a sequence "loses" the last two terms
        for m in range(len(A)-2):
            aitken_sequence.append(self.aitken(A, m))
        return (aitken_sequence)

    # Creates all transformed Aitken-sub-sequences in the class
    def set_aitken_normals(self):
        """
        Creates all transformed explicit Aitken-sub-sequences in the class
        
        """
        # Initial sub-sequence is the original sequence
        self.aitken_sequence_normal.append(self.aitken_normal(self.sequence))
        
        while (len(self.aitken_sequence_normal[-1]) >= 3):
            self.aitken_sequence_normal.append(self.aitken_normal(self.aitken_sequence_normal[-1]))

    # Aitken fixedpoint iteration
    def set_aitken_fixedpoint(self):
        """
        Creates a recursive sequence using the Steffensen/Aitken fixed-point iteration.
        Each iteration computes 2 function steps and uses the Aitken formula
        to reset the starting value for the next triplet.
        """
        # Initial value
        x_0 = self.sequence[0]
        
        # Maximum range of triplet steps
        steps = (self.max_m - self.min_m) // 3

        for _ in range(steps):
            # Next two fixed point steps
            x_1 = self.A(x_0)
            x_2 = self.A(x_1)
            
            # Save values
            self.aitken_sequence_fixedpoint.extend([x_0, x_1, x_2])
            
            # Denominator for Aitken
            denominator = x_2 - 2 * x_1 + x_0
            
            # Numeric stability test
            if abs(denominator) < 1e-15:
                break
                
            # Aitken-correction new initial value for next step
            x_0 = x_0 - ((x_1 - x_0) ** 2) / denominator

    # Create slides for presentation
    def show_slides(self, slides, bigdata=True):
        """
        Interactive presentation tool to visualize Aitken sequence transformations.

        Allows navigating through predefined slides using the arrow keys (left/right)
        and closing the viewer with the escape key. Plots the original sequence,
        the normally accelerated sub-sequences, and the fixed-point iteration.
        
        Input:
            - slides[list of tuples]: Definition of the presentation slides. Each tuple 
                                      contains: (title[str], show_lim[bool], show_og[bool], 
                                      normals[int], show_fixacc[bool])
            - bigdata[boolean]: If True, markers and legends are hidden to maintain 
                                performance and clarity for very large data sets.
        Output:
            - None: Displays an interactive Matplotlib figure window.
        """

        # State dictionary to keep track of the current slide index and presentation status.
        aktueller_index = {"i": 0, "running": True}

        # Initialize the matplotlib figure and axis for interactive plotting
        fig, ax = plt.subplots(figsize=(10, 6))
        fig.canvas.manager.set_window_title("Aitken Presentation")

        # Event handler for keyboard interactions
        def tastenevent(event):
            """Event handler for keyboard interactions."""
            # Advance to the next slide (wrap around to 0 if at the end)
            if event.key == "right":
                aktueller_index["i"] = (aktueller_index["i"] + 1) % len(slides)
                update_plot(bigdata)
            # Go back to the previous slide (wrap around to the last slide if at 0)
            elif event.key == "left":
                aktueller_index["i"] = (aktueller_index["i"] - 1) % len(slides)
                update_plot(bigdata)
            # Set running flag to False and close the current figure window
            elif event.key == "escape":
                aktueller_index["running"] = False
                plt.close(fig) 

        # Update slides
        def update_plot(bigdata):
            """Clears the previous frame and renders the current slide configuration."""

            # Clear the previous plot content
            ax.clear()  

            # Unpack the configurations for the current active slide
            idx = aktueller_index["i"]
            title, show_lim, show_og, normals, show_fixacc = slides[idx]

            # Adjust markers based on the data size to optimize performance and readability
            # No markers for large data sets
            if (bigdata):
                m_f = ""
                m_a = ""
            # Use visible markers for smaller data sets
            else:
                m_f = "o"
                m_a = "x"
            
            # Plot the hypothesized limit horizontal line if requested
            if show_lim and self.lim is not None:
                ax.axhline(self.lim, color='red', linestyle=':', linewidth=2, label=f"Hypothesized Limit: {self.lim}")
            
            # Plot the original sequence
            if show_og:
                m_range = list(range(self.min_m, self.max_m + 1))
                ax.plot(m_range, self.sequence, label="Original Sequence", linestyle='-', marker=m_f)
            
            # Plot normally accelerated Aitken sub-sequences based on the selection code
            # Option 0: Plot all generated sub-sequences
            if normals == 0:
                for i, folge in enumerate(self.aitken_sequence_normal):
                    m_aitken = list(range(self.min_m, self.min_m + len(folge)))
                    ax.plot(m_aitken, folge, label=f"Aitken normal {i+1} times", linestyle='--', marker=m_a)
            # Option 1: Plot only the first and the very last sub-sequence
            elif normals == 1:
                m_aitken = list(range(self.min_m, self.min_m + len(self.aitken_sequence_normal[0])))
                ax.plot(m_aitken, self.aitken_sequence_normal[0], label=f"Aitken normal 1 time", linestyle='--', marker=m_a)
                m_aitken = list(range(self.min_m, self.min_m + len(self.aitken_sequence_normal[-1])))
                ax.plot(m_aitken, self.aitken_sequence_normal[-1], label=f"Aitken normal {len(self.aitken_sequence_normal)} times", linestyle='--', marker=m_a)
            # Option 2: Plot only the very last sub-sequence
            elif normals == 2:
                m_aitken = list(range(self.min_m, self.min_m + len(self.aitken_sequence_normal[-1])))
                ax.plot(m_aitken, self.aitken_sequence_normal[-1], label=f"Aitken normal {len(self.aitken_sequence_normal)} times", linestyle='--', marker=m_a)
            # Option 3: Plot only the first sub-sequence
            elif normals == 3:
                m_aitken = list(range(self.min_m, self.min_m + len(self.aitken_sequence_normal[0])))
                ax.plot(m_aitken, self.aitken_sequence_normal[0], label=f"Aitken normal 1 time", linestyle='--', marker=m_a)
            
            # Plot the Steffensen/Aitken fixed-point accelerated sequence
            if show_fixacc:
                m_fixpunkt = list(range(self.min_m, self.min_m + len(self.aitken_sequence_fixedpoint)))
                ax.plot(m_fixpunkt, self.aitken_sequence_fixedpoint, label=f"Aitken Fixed-Point Accelerated", linestyle='--', marker=m_a)

            # Chart Styling & Labels
            ax.set_title(title, fontsize=14)
            ax.set_xlabel("Index m", fontsize=12)
            ax.set_ylabel("Sequence Terms", fontsize=12)
            ax.grid(True, linestyle=':', alpha=0.7)

            # Only display the legend if we are not dealing with massive datasets
            if (not bigdata):
                ax.legend()
            
            fig.tight_layout()

            # Request a redraw of the canvas without blocking
            fig.canvas.draw_idle()

        # Register the keypress event callback to the figure canvas
        fig.canvas.mpl_connect("key_press_event", tastenevent)

        # Render the initial slide upon opening
        update_plot(bigdata)

        # Display the window; execution blocks here until closed or ESC is pressed
        plt.show()

# Leibnitz series
def leibnitz_series(m):
    ln2 = 0
    for k in range(1, m+1):
        ln2 += (-1)**(k+1)/k
    return (ln2)


if (__name__ == "__main__"):
    # Global presentation settings
    extra = 0
    bigdata = False
    

    # =============================================================================== #
    # CASE 0: Sequence A_m = 1/m, Limit: 0.0 (Explicit Sequence)
    # =============================================================================== #
    folge_null = Sequence_Aitken(lambda m: 1/m, 0, 1, 12+extra, False, 1, False)

    # Compute the baseline sequence and apply normal Aitken acceleration loops
    folge_null.set_sequence_normal()
    folge_null.set_aitken_normals()


    # Define the slideshow presentation layout for this case
    name = "A_m = 1/m"
    slides = []
    # Slide 1: Original sequence and the limit line
    slides.append([name, True, True, -1, False])
    # Slide 2: Original sequence, limit line, and ALL generated normal Aitken sequences
    slides.append([name, True, True, 0, False])
    # Slide 3: Original sequence, limit line, plus only the first and last Aitken sequences
    slides.append([name, True, True, 1, False])

    # Console output and launch interactive viewer
    print("Null:")
    print(folge_null)
    folge_null.show_slides(slides, bigdata)
    # =============================================================================== #


    # =============================================================================== #
    # CASE 1: Calculating ln(2) using the Leibniz Series, Limit: ln(2)
    # This acts as an explicit sequence representing alternating partial sums.
    # =============================================================================== #
    folge_ln = Sequence_Aitken(leibnitz_series, log(2), 1, 8+extra, False, 1, False)

    # Compute the baseline sequence and apply normal Aitken acceleration loops
    folge_ln.set_sequence_normal()
    folge_ln.set_aitken_normals()

    # Define the slideshow presentation layout for this case
    name = "Leibnitzseries"
    slides = []
    # Slide 1: Original alternating sequence and the limit line
    slides.append([name, True, True, -1, False])
    # Slide 2: Original sequence, limit line, and ALL normal Aitken sequences
    slides.append([name, True, True, 0, False])
    # Slide 3: Original sequence, limit line, and only the FIRST Aitken sequence
    slides.append([name, True, True, 3, False])
    # Slide 4: Original sequence, limit line, and only the LAST Aitken sequence
    slides.append([name, True, True, 2, False])

    # Console output and launch interactive viewer
    print("ln:")
    print(folge_ln)
    folge_ln.show_slides(slides, bigdata)
    # =============================================================================== #


    # =============================================================================== #
    # CASE 2: Slow/Poorly behaving sequence: A_m = 1 - 0.70**m + 1/m, Limit: 1.0
    # =============================================================================== #
    series_bad = Sequence_Aitken(lambda m: 1 - 0.70**m + 1/m, 1, 1, 32+extra, False, 1, False)

    # Compute the baseline sequence and apply normal Aitken acceleration loops
    series_bad.set_sequence_normal()
    series_bad.set_aitken_normals()

    # Define the slideshow presentation layout for this case
    name = "A_m = 1 - 0.70**m + 1/m"
    slides = []
    # Slide 1: Original sequence and the limit line
    slides.append([name, True, True, -1, False])
    # Slide 2: Original sequence, limit line, and only the FIRST Aitken sequence
    slides.append([name, True, True, 3, False])
    # Slide 3: Original sequence, limit line, and only the LAST Aitken sequence
    slides.append([name, True, True, 2, False])

    # Console output and launch interactive viewer
    print("Bad:")
    print(series_bad)
    series_bad.show_slides(slides, bigdata)

    # =============================================================================== #


    # =============================================================================== #
    # CASE 3: Divergent / Oscillating Sequence: A_m = (-1)^m + 1/m
    # =============================================================================== #
    folge_divergent = Sequence_Aitken(lambda m: (-1)**m + 1/m, 0, 1, 16+extra, True, 1, False)

    # Compute the baseline sequence and apply normal Aitken acceleration loops
    folge_divergent.set_sequence_normal()
    folge_divergent.set_aitken_normals()

    # Define the slideshow presentation layout for this case
    name = "A_m = (-1)^m + 1/m"
    slides = []
    # Slide 1: Purely the original oscillating sequence
    slides.append([name, False, True, -1, False])
    # Slide 2: Original sequence along with ALL normal Aitken sequences
    slides.append([name, False, True, 0, False])
    # Slide 3: Original sequence along with only the FIRST Aitken sequence
    slides.append([name, False, True, 3, False])

    # Console output and launch interactive viewer
    print("Divergent:")
    print(folge_divergent)
    folge_divergent.show_slides(slides, bigdata)
    # =============================================================================== #


    # =============================================================================== #
    # CASE 4: Solving cos(x) = x via fixed-point iteration (Dottie Number)
    # Limit: ~0.739085 (Recursive Sequence)
    # =============================================================================== #
    folge_cos = Sequence_Aitken(cos, 0.7390851332151606, 0, 8+extra, True, 1, False)

    # Compute recursive sequence and apply both standard Aitken and Steffensen (Fixed-Point) loops
    folge_cos.set_sequence_recursive()
    folge_cos.set_aitken_normals()
    folge_cos.set_aitken_fixedpoint()

    # Define the slideshow presentation layout for this case
    name = "cos(x) = x"
    slides = []
    # Slide 1: Original recursive sequence and the target limit line
    slides.append([name, True, True, -1, False])
    # Slide 2: Original sequence, limit line, and ALL normal Aitken sequences
    slides.append([name, True, True, 0, False])
    # Slide 3: Original sequence, limit line, and only the LAST normal Aitken sequence
    slides.append([name, True, True, 2, False])
    # Slide 4: Original sequence, limit line, and the powerful Fixed-Point Acceleration curve
    slides.append([name, True, True, -1, True])
    # Slide 5: Clean comparison: Just the limit line, the last normal Aitken, and Fixed-Point Acceleration
    slides.append([name, True, False, 2, True])

    # Console output and launch interactive viewer
    print("Cos:")
    print(folge_cos)
    folge_cos.show_slides(slides, bigdata)
    # =============================================================================== #


    # =============================================================================== #
    # CASE 5: Solving arctan(x) = x via fixed-point iteration, Limit: 0.0
    # This case showcases extremely slow logarithmic convergence near 0.
    # =============================================================================== #
    folge_atan = Sequence_Aitken(lambda x: atan(x), 0, 0, 16+extra, True, 0.5, False)

    # Compute recursive sequence and apply both acceleration schemes
    folge_atan.set_sequence_recursive()
    folge_atan.set_aitken_normals()
    folge_atan.set_aitken_fixedpoint()

    # Define the slideshow presentation layout for this case
    name = "arctan(x) = x"
    slides = []
    # Slide 1: Original slow sequence and the limit line
    slides.append([name, True, True, -1, False])
    # Slide 2: Original sequence, limit line, and ALL normal Aitken sequences
    slides.append([name, True, True, 0, False])
    # Slide 3: Original sequence, limit line, and the rapid Fixed-Point Acceleration trajectory
    slides.append([name, True, True, -1, True])

    # Console output and launch interactive viewer
    print("Atan:")
    print(folge_atan)
    folge_atan.show_slides(slides, bigdata)
    # =============================================================================== #


    # =============================================================================== #
    # CASE 6: Solving cos(1/x) = x (Highly volatile/chaotic recursive behavior)
    # =============================================================================== #
    folge_chaos = Sequence_Aitken(lambda x: cos(1/x), 0, 0, 32+extra, True, 1, False)

    # Compute recursive sequence and apply both acceleration schemes
    folge_chaos.set_sequence_recursive()
    folge_chaos.set_aitken_normals()
    folge_chaos.set_aitken_fixedpoint()

    # Define the slideshow presentation layout for this case
    name = "cos(1/x) = x"
    slides = []
    # Slide 1: Purely the chaotic original sequence
    slides.append([name, False, True, -1, False])
    # Slide 2: Original sequence along with ALL normal Aitken lines
    slides.append([name, False, True, 0, False])
    # Slide 3: Original sequence along with only the FIRST Aitken sequence
    slides.append([name, False, True, 3, False])
    # Slide 4: Original sequence combined with the accelerated Fixed-Point path
    slides.append([name, False, True, -1, True])
    # Slide 5: Clean visualization of ONLY the accelerated Fixed-Point sequence tracking to the root
    slides.append([name, False, False, -1, True])

    # Print validation calculations to console to confirm correctness of the found root
    print("Chaos:") 
    print(folge_chaos)
    print("")
    print("Found root position x*:", folge_chaos.aitken_sequence_fixedpoint[-1])
    print("Verification cos(1/x*):", cos(1/folge_chaos.aitken_sequence_fixedpoint[-1]))
    
    # Launch interactive viewer for the final case
    folge_chaos.show_slides(slides, bigdata)
    # =============================================================================== #
