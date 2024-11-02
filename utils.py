def convert_to_float(input_str):
    try:
        # First try to convert directly to float (for inputs like "6" or "0")
        return float(input_str)
    except ValueError:
        # If that fails, try to handle it as a fraction
        if '/' in input_str:
            try:
                num, denom = input_str.split('/')
                return float(num) / float(denom)
            except (ValueError, ZeroDivisionError):
                raise ValueError("Invalid fraction format or division by zero")
        else:
            raise ValueError("Invalid number format")
