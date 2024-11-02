import tkinter as tk
from tkinter import ttk, messagebox
from ortools.linear_solver import pywraplp


class Optimiz(tk.Tk):
    def __init__(self):
        super().__init__()

        self.title("optimiz")
        self.geometry("600x500")

        # Problem Type section
        problem_frame = ttk.LabelFrame(self, text="Problem Type")
        problem_frame.pack(fill="x", padx=10, pady=5)
        ttk.Label(problem_frame, text="Problem Type").pack(side="left")
        self.problem_type = tk.StringVar(problem_frame)
        self.problem_type.set("")
        problem_type_options = ["-Choose-", "simple", "transportation", "hello"]
        ttk.OptionMenu(problem_frame, self.problem_type, *problem_type_options).pack(side="left")

        # Variables section
        self.var_frame = ttk.LabelFrame(self, text="Variables", padding="10")
        self.var_frame.pack(fill="x", padx=10, pady=5)

        ttk.Button(problem_frame, text="Set", command=self.setup_problem_scene).pack(side="left")

        # Objective function section
        obj_frame = ttk.LabelFrame(self, text="Objective Function", padding="10")
        obj_frame.pack(fill="x", padx=10, pady=5)

        self.obj_type = tk.StringVar(value="maximize")
        ttk.Radiobutton(obj_frame, text="Maximize", variable=self.obj_type,
                        value="maximize").pack(side="left")
        ttk.Radiobutton(obj_frame, text="Minimize", variable=self.obj_type,
                        value="minimize").pack(side="left")
        self.obj_coeffs_frame = ttk.Frame(obj_frame)
        self.obj_coeffs_frame.pack(fill="x", pady=5)

        # Constraints section
        const_frame = ttk.LabelFrame(self, text="Constraints", padding="10")
        const_frame.pack(fill="x", padx=10, pady=5)

        constraint_controls = ttk.Frame(const_frame)
        constraint_controls.pack(fill="x")

        ttk.Label(constraint_controls, text="Number of constraints:").pack(side="left")
        self.constraint_count = ttk.Entry(constraint_controls, width=5)
        self.constraint_count.pack(side="left", padx=5)
        ttk.Button(constraint_controls, text="Set",
                   command=self.setup_constraints).pack(side="left")

        self.constraints_frame = ttk.Frame(const_frame)
        self.constraints_frame.pack(fill="x", pady=5)

        # Solve button
        ttk.Button(self, text="Solve", command=self.solve).pack(pady=10)

        # Results section
        self.result_text = tk.Text(self, height=5, width=50)
        self.result_text.pack(pady=10, padx=10)


    def setup_problem_scene(self):
        if self.problem_type.get() == "simple":
            for widget in self.var_frame.winfo_children():
                widget.destroy()
            self.setup_simple_problem_scene()
        elif self.problem_type.get() == "transportation":
            for widget in self.var_frame.winfo_children():
                widget.destroy()
            self.setup_transportation_problem_scene()
        else:
            raise ValueError("We're under development")

    def setup_simple_problem_scene(self):
        ttk.Label(self.var_frame, text="Number of variables:").pack(side="left")
        self.var_count = ttk.Entry(self.var_frame, width=5)
        self.var_count.pack(side="left", padx=5)
        ttk.Button(self.var_frame, text="Set", command=self.setup_simple_problem_var).pack(side="left")


    def setup_transportation_problem_scene(self):
        origin_frame = (ttk.LabelFrame(self.var_frame, text="Origins"))
        origin_frame.grid(row=0, column=0)
        # Set up headers
        tk.Label(origin_frame, text="Origin").grid(row=0, column=0, padx=5, pady=5)
        tk.Label(origin_frame, text="Supply").grid(row=0, column=1, padx=5, pady=5)
        tk.Label(origin_frame, text="Cost").grid(row=0, column=2, padx=5, pady=5)
        for i in range(3):
            tk.Label(origin_frame, text=f"Origin {i+1}").grid(row=i+1, column=0, padx=5, pady=5)


        dest_frame = ttk.LabelFrame(self.var_frame, text="Destinations")
        dest_frame.grid(row=0, column=1)
        relative_col_idx = 3
        tk.Label(dest_frame, text="Destination").grid(row=0, column=relative_col_idx, padx=5, pady=5)
        tk.Label(dest_frame, text="Demand").grid(row=0, column=relative_col_idx+1, padx=5, pady=5)

        for i in range(3):
            tk.Label(dest_frame, text=f"Destination {i+1}").grid(row=i+1, column=relative_col_idx, padx=5, pady=5)

    #label.pack(side="left")
        # ttk.Label(self.var_frame, text="Number of variables:").pack(side="left")
        # self.var_count = ttk.Entry(self.var_frame, width=5)
        # self.var_count.pack(side="left", padx=5)
        # ttk.Button(self.var_frame, text="Set", command=self.setup_simple_problem_var).pack(side="left")

    def setup_simple_problem_var(self):
        try:
            n = int(self.var_count.get())
            if n <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Error", "Please enter a positive integer")
            return

        for widget in self.obj_coeffs_frame.winfo_children():
            widget.destroy()

        for i in range(n):
            ttk.Label(self.obj_coeffs_frame, text=f"x{i+1}:").pack(side="left")
            entry = ttk.Entry(self.obj_coeffs_frame, width=5)
            entry.pack(side="left", padx=2)

    def setup_constraints(self):
        try:
            n = int(self.constraint_count.get())
            var_count = len(self.obj_coeffs_frame.winfo_children()) // 2
            if n <= 0 or var_count == 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Error", "Please set valid number of variables and constraints")
            return

        for widget in self.constraints_frame.winfo_children():
            widget.destroy()

        for i in range(n):
            constraint_frame = ttk.Frame(self.constraints_frame)
            constraint_frame.pack(fill="x", pady=2)

            for j in range(var_count):
                ttk.Label(constraint_frame, text=f"x{j+1}:").pack(side="left")
                entry = ttk.Entry(constraint_frame, width=5)
                entry.pack(side="left", padx=2)

            relation = ttk.Combobox(constraint_frame, values=["<=", ">=", "="], width=3)
            relation.set("<=")
            relation.pack(side="left", padx=5)

            rhs = ttk.Entry(constraint_frame, width=5)
            rhs.pack(side="left", padx=2)

    def solve(self):
        try:
            # Create solver
            solver = pywraplp.Solver.CreateSolver('GLOP')
            if not solver:
                messagebox.showerror("Error", "Could not create solver")
                return

            # Get number of variables
            var_count = len(self.obj_coeffs_frame.winfo_children()) // 2

            # Create variables
            variables = []
            for i in range(var_count):
                variables.append(solver.NumVar(0, solver.infinity(), f'x{i+1}'))

            # Set objective function
            obj_coeffs = []
            for i in range(0, len(self.obj_coeffs_frame.winfo_children()), 2):
                entry = self.obj_coeffs_frame.winfo_children()[i+1]
                obj_coeffs.append(float(entry.get()))

            objective = solver.Objective()
            for i, coeff in enumerate(obj_coeffs):
                objective.SetCoefficient(variables[i], coeff)

            if self.obj_type.get() == "maximize":
                objective.SetMaximization()
            else:
                objective.SetMinimization()

            # Add constraints
            for constraint_frame in self.constraints_frame.winfo_children():
                widgets = constraint_frame.winfo_children()
                coeffs = []
                for i in range(0, var_count * 2, 2):
                    entry = widgets[i+1]
                    coeffs.append(float(entry.get()))

                relation = widgets[-2].get()
                rhs = float(widgets[-1].get())

                constraint = solver.Constraint(-solver.infinity(), solver.infinity())
                for i, coeff in enumerate(coeffs):
                    constraint.SetCoefficient(variables[i], coeff)

                if relation == "<=":
                    constraint.SetUB(rhs)
                elif relation == ">=":
                    constraint.SetLB(rhs)
                else:  # "="
                    constraint.SetBounds(rhs, rhs)

            # Solve
            status = solver.Solve()

            # Display results
            self.result_text.delete(1.0, tk.END)
            if status == pywraplp.Solver.OPTIMAL:
                self.result_text.insert(tk.END, f"Solution found!\n")
                self.result_text.insert(tk.END, f"Objective value = {solver.Objective().Value()}\n")
                for i, var in enumerate(variables):
                    self.result_text.insert(tk.END, f"x{i+1} = {var.solution_value()}\n")
            else:
                self.result_text.insert(tk.END, "The problem does not have an optimal solution.")

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")


app = Optimiz()
app.mainloop()