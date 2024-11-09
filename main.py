import tkinter as tk
from tkinter import ttk, messagebox
from ortools.linear_solver import pywraplp
from utils import convert_to_float


class SimpleFrame(ttk.LabelFrame):
    def __init__(self, master):
        super().__init__(master)

        # Variables section
        var_frame = ttk.LabelFrame(self, text="Variables", padding="10")
        var_frame.pack(fill="x", padx=10, pady=5)
        self.n_var = tk.IntVar(value=1)
        self.setup_meta_variables(var_frame, self.n_var)
        # ttk.Button(var_frame, text="Set", command=self.setup_variables).pack(side="left")

        # Objective function section
        objective_frame = ttk.LabelFrame(self, text="Objective Function", padding="10")
        objective_frame.pack(fill="x", padx=10, pady=5)

        # Objective type
        self.objective_type = tk.StringVar(value="maximize")
        ttk.Radiobutton(objective_frame, text="Maximize", variable=self.objective_type,
                        value="maximize").pack(side="left")
        ttk.Radiobutton(objective_frame, text="Minimize", variable=self.objective_type,
                        value="minimize").pack(side="left")
        self.objective_coeffs_frame = ttk.Frame(objective_frame)
        self.objective_coeffs_frame.pack(side="bottom", fill="x", pady=5)
        self.objective_coeffs = []

        # Constraints section
        const_frame = ttk.LabelFrame(self, text="Constraints", padding="10")
        const_frame.pack(fill="x", padx=10, pady=5)
        cons_controls_frame = ttk.Frame(const_frame)
        cons_controls_frame.pack(fill="x")
        ttk.Label(cons_controls_frame, text="Number of constraints:").pack(side="left")
        self.n_constraints = tk.IntVar(value=1)
        ttk.Entry(cons_controls_frame, width=5, textvariable=self.n_constraints).pack(side="left", padx=5)
        ttk.Button(cons_controls_frame, text="Set", command=self.setup_constraints_equations).pack(side="left")
        self.const_eq_frame = ttk.Frame(const_frame)
        self.const_eq_frame.pack(fill="x", pady=5)
        self.const_names = [] # Optional name for each coeffs
        self.const_coeffs = []  # const_coeffs[const_num][coeff_num]  const_num in [0, n_coeffs), coeff_num in [0, n_var)
        self.const_relations = [] # one relation per constraint equation "<=, ==, >="
        self.const_rhs = []

        # Solve button
        ttk.Button(self, text="Solve", command=self.solve).pack(pady=10)

        # Results section
        self.result_text = tk.Text(self, height=25, width=75)
        self.result_text.pack(pady=10, padx=10)

    def setup_meta_variables(self, var_frame, variable_count):
        ttk.Label(var_frame, text="Number of variables:").pack(side="left")
        ttk.Entry(var_frame, width=5, textvariable=variable_count).pack(side="left", padx=5)
        ttk.Button(var_frame, text="Set", command=self.setup_objective_coeffs).pack(side="left", padx=5)

    def setup_objective_coeffs(self):
        try:
            if self.n_var.get() <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Error", "Please enter a positive integer")
            return

        for widget in self.objective_coeffs_frame.winfo_children():
            widget.destroy()
        self.objective_coeffs = []

        for i in range(self.n_var.get()):
            coeff_i = tk.IntVar(value=0)
            ttk.Label(self.objective_coeffs_frame, text=f"x{i+1}:").pack(side="left")
            ttk.Entry(self.objective_coeffs_frame, width=5, textvariable=coeff_i).pack(side="left", padx=2)
            self.objective_coeffs.append(coeff_i)

    def setup_constraints_equations(self):
        try:
            if self.n_constraints.get() <= 0 or self.n_var.get() <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Error", "Please set valid number of variables and constraints")
            return

        for widget in self.const_eq_frame.winfo_children():
            widget.destroy()
        self.const_names = []
        self.const_coeffs = []
        self.const_relations = []
        self.const_rhs = []

        for i in range(self.n_constraints.get()):
            coeffs_eq_i = []
            eq_frame = ttk.Frame(self.const_eq_frame)
            eq_frame.pack(fill="x", pady=2)
            const_name = tk.StringVar(value=f"const_{i+1}")
            ttk.Entry(eq_frame, width=10, name=f"const_{i+1}_name", textvariable=const_name).pack(side="left", padx=5)
            self.const_names.append(const_name)

            for j in range(self.n_var.get()):
                ttk.Label(eq_frame, text=f"x{j+1}:").pack(side="left")
                coeff_i_j = tk.StringVar(value="0")
                ttk.Entry(eq_frame, width=5, name=f"constraint_{i+1}_{j+1}_coeff", textvariable=coeff_i_j).pack(side="left", padx=2)
                coeffs_eq_i.append(coeff_i_j)
            self.const_coeffs.append(coeffs_eq_i)

            relation = tk.StringVar(value="<=")
            ttk.Combobox(eq_frame, values=["<=", ">=", "="], width=3, textvariable=relation).pack(side="left", padx=2)
            self.const_relations.append(relation)

            rhs = tk.StringVar(value="0")
            ttk.Entry(eq_frame, width=5, textvariable=rhs).pack(side="left", padx=2)
            self.const_rhs.append(rhs)

    def solve(self):
        constraints_info_dict = dict() # name: {expression: '3x+4y<=6', solution_value: 1337, slack/surplus: 0}
        # Create solver
        solver = pywraplp.Solver.CreateSolver('GLOP')
        if not solver:
            messagebox.showerror("Error", "Could not create solver")
            return

        # Create variables
        variables = []
        for i in range(self.n_var.get()):
            variables.append(solver.NumVar(0, solver.infinity(), f'x{i+1}'))

        # Set objective function
        objective = solver.Objective()
        for i, coeff in enumerate(self.objective_coeffs):
            objective.SetCoefficient(variables[i], convert_to_float(coeff.get()))

        if self.objective_type.get() == "maximize":
            objective.SetMaximization()
        else:
            objective.SetMinimization()

        # Add constraints
        constraints = []
        for const_idx in range(self.n_constraints.get()):
            constraint_name = self.const_names[const_idx].get()
            constraints_info_dict[constraint_name] = dict()
            coeffs = []
            coeffs_str = []
            for var_idx in range(self.n_var.get()):
                coeff = self.const_coeffs[const_idx][var_idx].get()
                coeffs.append(convert_to_float(coeff))
                coeffs_str.append(coeff)

            relation = self.const_relations[const_idx].get()
            rhs_str = self.const_rhs[const_idx].get()
            rhs = convert_to_float(rhs_str)

            constraints_info_dict[constraint_name]["expr"] = " + ".join(coeffs_str) + f" {relation} {rhs_str}"
            print(constraints_info_dict[constraint_name]["expr"])
            constraint = solver.Constraint(-solver.infinity(), solver.infinity())
            for i, coeff in enumerate(coeffs):
                constraint.SetCoefficient(variables[i], coeff)

            if relation == "<=":
                constraint.SetUb(rhs)
                constraints_info_dict[constraint_name]["slack"] = "tofill"
            elif relation == ">=":
                constraint.SetLb(rhs)
                constraints_info_dict[constraint_name]["surplus"] = "tofill"
            else:  # "="
                constraint.SetBounds(rhs, rhs)
            constraints.append(constraint)

        # Solve
        status = solver.Solve()

        # Display results
        self.result_text.delete(1.0, tk.END)
        if status == pywraplp.Solver.OPTIMAL:
            self.result_text.insert(tk.END, f"Solution found!\n\n")
            self.result_text.insert(tk.END, f"Objective value = {solver.Objective().Value()}\n")
            for i, var in enumerate(variables):
                self.result_text.insert(tk.END, f"x{i+1} = {var.solution_value()}\n")

            self.result_text.insert(tk.END, f"\n ----------------- \n\n")

            def compute_slack_or_surplus(constraint):
                constraint_solution_value = 0
                for v in variables:
                    constraint_solution_value += v.solution_value() * constraint.GetCoefficient(v)
                if constraint.ub() != solver.infinity():
                    # the upper bound is set so it will give slack if any
                    slack = constraint.ub() - constraint_solution_value
                    return constraint_solution_value, slack
                else:
                    surplus = constraint_solution_value - constraint.lb()
                    return constraint_solution_value, surplus

            for constraint_name, constraint in zip(constraints_info_dict, constraints):
                const_str = ""
                # const_str += constraint_name + ": " + f"{constraints_info_dict[constraint_name]['expr']}"
                const_str += constraint_name + ": "
                if "slack" in constraints_info_dict[constraint_name]:
                    solution_value, slack = compute_slack_or_surplus(constraint)
                    const_str += f" solution={solution_value}  slack={slack}"
                elif "surplus" in constraints_info_dict[constraint_name]:
                    solution_value, surplus = compute_slack_or_surplus(constraint)
                    const_str += f" solution={solution_value}  surplus={surplus}"
                else:
                    const_str += f" slack/surplus = 0"
                self.result_text.insert(tk.END, const_str + "\n")
        else:
            self.result_text.insert(tk.END, "The problem does not have an optimal solution.")


class Optimiz(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("optimiz")
        self.geometry("1200x1000")

        # Problem Type section
        problem_type_frame = ttk.LabelFrame(self, text="Problem Type")
        problem_type_frame.pack(fill="x", padx=10, pady=1)
        ttk.Label(problem_type_frame, text="Problem Type").pack(side="left")
        self.problem_type = tk.StringVar(problem_type_frame)
        self.problem_type.set("")

        problem_type_options = ["--Choose--", "simple", "transportation", "hello"]
        ttk.OptionMenu(problem_type_frame, self.problem_type, *problem_type_options).pack(side="left")
        ttk.Button(problem_type_frame, text="Set", command=self.setup_problem_scene).pack(side="left")

        self.problem_frame = None

    def setup_problem_scene(self):
        if self.problem_type.get() == "simple":
            self.problem_frame = SimpleFrame(self)
            self.problem_frame.pack(fill="both", side="left", expand=True)

    def setup_transportation_problem_scene(self):
        ttk.Label(self.var_frame, text="Num Origin", name="n_origins").grid(row=0, column=0)
        self.n_origins = ttk.Entry(self.var_frame, width=5)
        self.n_origins.grid(row=0, column=1, padx=5)
        ttk.Label(self.var_frame, text="Num Destinations", name="n_destinations").grid(row=0, column=2)
        self.n_destinations = ttk.Entry(self.var_frame, width=5)
        self.n_destinations.grid(row=0, column=3, padx=5)
        ttk.Button(self.var_frame, text="Set", command=self.setup_variables_transportation).grid(row=0, column=4)

    #label.pack(side="left")
        # ttk.Label(self.var_frame, text="Number of variables:").pack(side="left")
        # self.n_var = ttk.Entry(self.var_frame, width=5)
        # self.n_var.pack(side="left", padx=5)
        # ttk.Button(self.var_frame, text="Set", command=self.setup_simple_problem_var).pack(side="left")

    def setup_variables_transportation(self):
        try:
            n_origins = int(self.n_origins.get())
            n_destinations = int(self.n_destinations.get())
            if n_origins <= 0 or n_destinations <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Error", "Please enter a positive integer")
            return

        self.transportation_vars_frame = ttk.Frame(self.var_frame)
        self.transportation_vars_frame.grid(row=1)

        # Setup headers
        for i in range(n_origins):
            ttk.Label(self.transportation_vars_frame, text=f"Origin {i}").grid(row=i+1, column=0)
        for j in range(n_destinations):
            ttk.Label(self.transportation_vars_frame, text=f"Destination {j}").grid(row=0, column=j+1)


        #for row in range(n_destinations + 1):
        # origin_frame = ttk.LabelFrame(self.var_frame, text="Origins")
        # tk.Label(origin_frame, text="Supply").grid(row=0, column=1, padx=5, pady=5)
        # tk.Label(origin_frame, text="Cost").grid(row=0, column=2, padx=5, pady=5)
        # for i in range(3):
        #     tk.Label(origin_frame, text=f"Origin {i+1}").grid(row=i+1, column=0, padx=5, pady=5)
        #
        # dest_frame = ttk.LabelFrame(self.var_frame, text="Destinations")
        # dest_frame.grid(row=0, column=1)
        # relative_col_idx = 3
        # tk.Label(dest_frame, text="Destination").grid(row=0, column=relative_col_idx, padx=5, pady=5)
        # tk.Label(dest_frame, text="Demand").grid(row=0, column=relative_col_idx+1, padx=5, pady=5)
        #
        # for i in range(3):
        #     tk.Label(dest_frame, text=f"Destination {i+1}").grid(row=i+1, column=relative_col_idx, padx=5, pady=5)



app = Optimiz()
app.mainloop()