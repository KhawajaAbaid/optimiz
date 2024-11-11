import tkinter as tk
from tkinter import ttk, messagebox
from ortools.linear_solver import pywraplp
from utils import convert_to_float
from functools import partial


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
            # print(constraints_info_dict[constraint_name]["expr"])
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


class TransportationFrame(ttk.LabelFrame):
    def __init__(self, master):
        super().__init__(master)
        # Variables to store node counts
        self.n_sources = tk.IntVar(value=2)
        self.n_destinations = tk.IntVar(value=2)
        self.objective_type = tk.StringVar(value="minimize")

        # Create main frames
        self.setup_network_definition()
        self.setup_network_details()
        self.setup_controls()

        self.src_names = []
        self.src_supplies = []
        self.dest_names = []
        self.dest_demands = []
        self.cost_matrix = []   # cost_matrix[source_i][dest_j]

        # Results section
        self.result_text = tk.Text(self, height=25, width=75)
        self.result_text.pack(pady=10, padx=10)


    def clear_data_lists(self):
        self.src_names = []
        self.src_supplies = []
        self.dest_names = []
        self.dest_demands = []
        self.cost_matrix = []

    def setup_network_definition(self):
        """Create frame for node count inputs"""
        node_frame = ttk.LabelFrame(self, text="Network Definition", padding="10")
        node_frame.pack(fill="x", padx=10, pady=5)

        # Source nodes
        ttk.Label(node_frame, text="Number of Sources:").grid(row=0, column=0, padx=5, pady=5)
        ttk.Entry(node_frame, textvariable=self.n_sources, width=10).grid(row=0, column=1, padx=5, pady=5)

        # Destination nodes
        ttk.Label(node_frame, text="Number of Destinations:").grid(row=0, column=4, padx=5, pady=5)
        ttk.Entry(node_frame, textvariable=self.n_destinations, width=10).grid(row=0, column=5, padx=5, pady=5)
        ttk.Button(node_frame, text="Set", command=self.update_tables).grid(row=0, column=6, padx=5, pady=5)
        self.clear_data_lists()

    def setup_network_details(self):
        """Create frame for cost inputs"""
        self.data_frame = ttk.Notebook(self)
        self.data_frame.pack(fill="both", expand=True, padx=10, pady=5)

        # Create tabs for different input types
        self.supply_frame = ttk.Frame(self.data_frame)
        self.demand_frame = ttk.Frame(self.data_frame)
        self.costs_frame = ttk.Frame(self.data_frame)
        self.results_frame = ttk.Frame(self.data_frame)

        self.data_frame.add(self.supply_frame, text="Supply")
        self.data_frame.add(self.demand_frame, text="Demand")
        self.data_frame.add(self.costs_frame, text="Transportation Costs")
        self.data_frame.add(self.results_frame, text="Results")

    def setup_controls(self):
        """Create frame for action buttons"""
        # Let's add objective before buttons
        objective_frame = ttk.Frame(self)
        objective_frame.pack(fill="x", padx=10, pady=5)
        ttk.Radiobutton(objective_frame, text="Maximize", variable=self.objective_type,
                        value="maximize").pack(side="left")
        ttk.Radiobutton(objective_frame, text="Minimize", variable=self.objective_type,
                        value="minimize").pack(side="left")

        button_frame = ttk.Frame(self)
        button_frame.pack(fill="x", padx=10, pady=5)
        ttk.Button(button_frame, text="Update Tables", command=self.update_tables).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Solve", command=self.solve).pack(side="left", padx=5)
        ttk.Button(button_frame, text="Clear", command=self.clear).pack(side="left", padx=5)

    def setup_costs(self):
        self.cost_matrix = []
        ttk.Label(self.costs_frame, text="Transportation Costs").pack(pady=5)
        for widget in self.costs_frame.winfo_children():
            widget.destroy()
        # Create cost matrix
        cost_matrix_frame = ttk.Frame(self.costs_frame)
        cost_matrix_frame.pack(fill="x", padx=5, pady=2)
        # -> setup destination headers
        for j in range(self.n_destinations.get()):
            dest_label = ttk.Label(cost_matrix_frame, text=self.dest_names[j].get())
            dest_label.grid(row=0, column=j+1, padx=5, pady=2)

        for i in range(self.n_sources.get()):
            src_label = ttk.Label(cost_matrix_frame, text=self.src_names[i].get())
            src_label.grid(row=i+1, column=0, padx=5, pady=2)
            src_i_costs = []
            for j in range(self.n_destinations.get()):
                cost_i_j = tk.StringVar(value="0")
                ttk.Entry(cost_matrix_frame, width=8, textvariable=cost_i_j).grid(row=i+1, column=j+1, padx=5, pady=2)
                src_i_costs.append(cost_i_j)
            self.cost_matrix.append(src_i_costs)

    def display_results(self, solution_matrix):
        for widget in self.results_frame.winfo_children():
            widget.destroy()
        cost_matrix_float = []
        for i in range(self.n_sources.get()):
            row = []
            for j in range(self.n_destinations.get()):
                row.append(convert_to_float(self.cost_matrix[i][j].get()))
            cost_matrix_float.append(row)

        # Create units shipped matrix
        ttk.Label(self.results_frame, text="Units Shipped").pack(pady=5)
        units_shipped_frame = ttk.Frame(self.results_frame)
        units_shipped_frame.pack(fill="x", padx=5, pady=2)
        # -> setup destination headers
        for j in range(self.n_destinations.get()):
            dest_label = ttk.Label(units_shipped_frame, text=self.dest_names[j].get())
            dest_label.grid(row=0, column=j+1, padx=5, pady=2)

        for i in range(self.n_sources.get()):
            src_label = ttk.Label(units_shipped_frame, text=self.src_names[i].get())
            src_label.grid(row=i+1, column=0, padx=5, pady=2)
            for j in range(self.n_destinations.get()):
                val_label = ttk.Label(units_shipped_frame, text=f"{solution_matrix[i][j]}")
                val_label.grid(row=i+1, column=j+1, padx=5, pady=2)

        ttk.Label(self.results_frame, text="Total Cost Per Route").pack(pady=5)
        total_cost = 0
        total_cost_per_route_frame = ttk.Frame(self.results_frame)
        total_cost_per_route_frame.pack(fill="x", padx=5, pady=2)
        # -> setup destination headers
        for j in range(self.n_destinations.get()):
            dest_label = ttk.Label(total_cost_per_route_frame, text=self.dest_names[j].get())
            dest_label.grid(row=0, column=j+1, padx=5, pady=2)

        for i in range(self.n_sources.get()):
            src_label = ttk.Label(total_cost_per_route_frame, text=self.src_names[i].get())
            src_label.grid(row=i+1, column=0, padx=5, pady=2)
            for j in range(self.n_destinations.get()):
                route_total_cost = solution_matrix[i][j] * cost_matrix_float[i][j]
                total_cost += route_total_cost
                val_label = ttk.Label(total_cost_per_route_frame, text=f"{route_total_cost}")
                val_label.grid(row=i+1, column=j+1, padx=5, pady=2)

        ttk.Label(self.results_frame, text=f"Total Cost: {total_cost}").pack(pady=5)

    def update_tables(self):
        """Update input tables based on node counts"""
        self.clear_data_lists()

        # Create supply inputs
        ttk.Label(self.supply_frame, text="Source Supply Quantities").pack(pady=5)
        for i in range(self.n_sources.get()):
            frame = ttk.Frame(self.supply_frame)
            frame.pack(fill="x", padx=5, pady=2)
            name = tk.StringVar(value=f"src_{i+1}")
            ttk.Entry(frame, textvariable=name).pack(side="left", padx=5)
            self.src_names.append(name)
            supply = tk.IntVar(value=0)
            ttk.Entry(frame, width=10, textvariable=supply).pack(side="left", padx=5)
            self.src_supplies.append(supply)
        ttk.Button(self.supply_frame, text="Set", command=self.setup_costs).pack(padx=5)

        # Create demand inputs
        ttk.Label(self.demand_frame, text="Destination Demand Quantities").pack(pady=5)
        for i in range(self.n_destinations.get()):
            frame = ttk.Frame(self.demand_frame)
            frame.pack(fill="x", padx=5, pady=2)
            name = tk.StringVar(value=f"dest_{i+1}")
            ttk.Entry(frame, textvariable=name).pack(side="left", padx=5)
            self.dest_names.append(name)
            demand = tk.IntVar(value=0)
            ttk.Entry(frame, width=10, textvariable=demand).pack(side="left", padx=5)
            self.dest_demands.append(demand)
        ttk.Button(self.demand_frame, text="Set", command=self.setup_costs).pack(padx=5)

        self.setup_costs()

    def solve(self):
        """Collect inputs and solve the transporation problem"""
        constraints_info_dict = dict() # name: {expression: '3x+4y<=6', solution_value: 1337, slack/surplus: 0}
        # Create solver
        solver = pywraplp.Solver.CreateSolver('GLOP')
        if not solver:
            messagebox.showerror("Error", "Could not create solver")
            return

        # Create variables
        var_matrix = []
        for i in range(self.n_sources.get()):
            src_i_var = []
            for j in range(self.n_destinations.get()):
                src_i_var.append(solver.NumVar(0, solver.infinity(), f'x{i+1}{j+1}'))
            var_matrix.append(src_i_var)

        # Set objective function
        # The values in cost matrix act as coeffs for objecive function
        objective = solver.Objective()
        for i in range(self.n_sources.get()):
            for j in range(self.n_destinations.get()):
                coeff = convert_to_float(self.cost_matrix[i][j].get())
                objective.SetCoefficient(var_matrix[i][j], coeff)

        if self.objective_type.get() == "minimize":
            objective.SetMinimization() # Do we need max?, well turns out we do
        else:
            objective.SetMaximization()

        # Set constraints
        constraints = []
        # Add supply constraints
        # The sum of quantities originating from each source should be less than or equal to the supply.
        # i.e. For a given source i, Sum[x_i_j for j in range(n_destinations)] <= supplies[i]
        for i in range(self.n_sources.get()):
            constraint_name = self.src_names[i].get() + '_quantity'
            constraint = solver.Constraint(-solver.infinity(), solver.infinity(), f'{constraint_name}_const')
            constraints_info_dict[constraint_name] = dict()
            coeffs_str = []
            for j in range(self.n_destinations.get()):
                coeff = 1
                constraint.SetCoefficient(var_matrix[i][j], 1)
                coeffs_str.append(str(coeff) + f"x{i+1}{j+1}")

            constraints_info_dict[constraint_name]["expr"] = " + ".join(coeffs_str) + f" <= {self.src_supplies[i].get()}"
            constraint.SetUb(convert_to_float(self.src_supplies[i].get()))
            constraints_info_dict[constraint_name]["slack"] = "tofill"
            constraints.append(constraint)

        # Add demand constraints
        # The sum of quantities received by the destination should be equal to the demand
        # i.e. For a given destination j, Sum[x_i_j for i in range(n_destinations)] = demands[j]
        for j in range(self.n_destinations.get()):
            constraint_name = self.dest_names[j].get() + '_quantity'
            constraint = solver.Constraint(-solver.infinity(), solver.infinity(), f'{constraint_name}_const')
            constraints_info_dict[constraint_name] = dict()
            coeffs_str = []
            for i in range(self.n_sources.get()):
                coeff = 1
                constraint.SetCoefficient(var_matrix[i][j], coeff)
                coeffs_str.append(str(coeff) + f"x{i+1}{j+1}")
            constraints_info_dict[constraint_name]["expr"] = " + ".join(coeffs_str) + f" = {self.dest_demands[j].get()}"
            demand = convert_to_float(self.dest_demands[j].get())
            constraint.SetBounds(demand, demand)
            constraints.append(constraint)

        # Solve
        status = solver.Solve()

        # Display results
        self.result_text.delete(1.0, tk.END)
        if status == pywraplp.Solver.OPTIMAL:
            self.result_text.insert(tk.END, f"Solution found!\n\n")
            self.result_text.insert(tk.END, f"Objective value = {solver.Objective().Value()}\n")
            solution_matrix = []
            for i in range(self.n_sources.get()):
                row = []
                for j in range(self.n_destinations.get()):
                    solution = round(var_matrix[i][j].solution_value(), 2)
                    row.append(solution)
                    self.result_text.insert(tk.END, f"x{i+1}{j} = {solution} \t")
                solution_matrix.append(row)
                self.result_text.insert(tk.END, "\n")
            self.result_text.insert(tk.END, f"\n ----------------- \n\n")
            self.display_results(solution_matrix)
        else:
            self.result_text.insert(tk.END, "The problem does not have an optimal solution.")

    def clear(self):
        """Clear all inputs"""
        self.n_sources.set(2)
        self.n_destinations.set(2)
        # Clear existing tables
        for widget in self.supply_frame.winfo_children():
            widget.destroy()
        for widget in self.demand_frame.winfo_children():
            widget.destroy()
        for widget in self.costs_frame.winfo_children():
            widget.destroy()
        self.update_tables()


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
        if self.problem_frame is not None:
            self.problem_frame.destroy()
        if self.problem_type.get() == "simple":
            self.problem_frame = SimpleFrame(self)
        if self.problem_type.get() == "transportation":
            self.problem_frame = TransportationFrame(self)

        self.problem_frame.pack(fill="both", side="left", expand=True)


    # def setup_transportation_problem_scene(self):
    #     ttk.Label(self.var_frame, text="Num Origin", name="n_origins").grid(row=0, column=0)
    #     self.n_origins = ttk.Entry(self.var_frame, width=5)
    #     self.n_origins.grid(row=0, column=1, padx=5)
    #     ttk.Label(self.var_frame, text="Num Destinations", name="n_destinations").grid(row=0, column=2)
    #     self.n_destinations = ttk.Entry(self.var_frame, width=5)
    #     self.n_destinations.grid(row=0, column=3, padx=5)
    #     ttk.Button(self.var_frame, text="Set", command=self.setup_variables_transportation).grid(row=0, column=4)

    #label.pack(side="left")
        # ttk.Label(self.var_frame, text="Number of variables:").pack(side="left")
        # self.n_var = ttk.Entry(self.var_frame, width=5)
        # self.n_var.pack(side="left", padx=5)
        # ttk.Button(self.var_frame, text="Set", command=self.setup_simple_problem_var).pack(side="left")

    # def setup_variables_transportation(self):
    #     try:
    #         n_origins = int(self.n_origins.get())
    #         n_destinations = int(self.n_destinations.get())
    #         if n_origins <= 0 or n_destinations <= 0:
    #             raise ValueError
    #     except ValueError:
    #         messagebox.showerror("Error", "Please enter a positive integer")
    #         return
    #
    #     self.transportation_vars_frame = ttk.Frame(self.var_frame)
    #     self.transportation_vars_frame.grid(row=1)
    #
    #     # Setup headers
    #     for i in range(n_origins):
    #         ttk.Label(self.transportation_vars_frame, text=f"Origin {i}").grid(row=i+1, column=0)
    #     for j in range(n_destinations):
    #         ttk.Label(self.transportation_vars_frame, text=f"Destination {j}").grid(row=0, column=j+1)


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