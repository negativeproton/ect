import sys

import curves


def main(print_welcome):
    if print_welcome:
        # Welcome message
        print("Welcome. This program transforms an elliptic curve of one type (and given coordinates) in another type.")
        print("Below you find the keys for the available types.")

    # Define an individual key for each curve type.
    curves_types = {
        'sw': curves.ShortWeierstrass,
        'm': curves.Montgomery,
        'ed': curves.Edwards,
        'ted': curves.TwistedEdwards
    }

    # Output all keys and the corresponding type.
    for key in curves_types:
        print("'" + key + "'" + "\t for " + curves_types[key].curve_name)

    # Includes input validation.
    user_input = get_user_input_base_key(curves_types)

    # Now, after input validation, it should be guaranteed, that user_input is a valid, specific key for a curve.

    if user_input == "ed":
        print("Currently there are only transformations that require c = 1. Please choose 1 for c.")

    # Dic for necessary arguments of all curves.
    # Values will get overwritten with user inputs.
    curve_arguments_dic = {
        "p": None,
        curves_types[user_input].coefficient1_name: None,
        curves_types[user_input].coefficient2_name: None,
        curves_types[user_input].coordinate1_name: None,
        curves_types[user_input].coordinate2_name: None
    }

    # Idea for future:
    # Output the possibility of default values for x and y when implemented here or in the following for loop with if statement key == "x".
    # Can be done simply by writing x=value in __init__. 

    for key in curve_arguments_dic.keys():
        # key should be a string.
        assert type(key) == str
        # Idea: advanced input validation with regex, only numbers, limit input size, avoid big numbers as power.
        # Allows to enter potency.

        # without error_handling function:
        # curve_arguments_dic[key] = eval(input("Please enter the value for " + key + ": "))

        # Calls input().
        curve_arguments_dic[key] = eval_parameter_input(key)

        assert type(curve_arguments_dic[key]) == int

    # The values should now be written in the curve_arguments_dic to the right key. 

    user_curve = init_curve(curve_arguments_dic, curves_types, user_input)

    # Now there should be the base curve initialized called user_curve.

    ask_transtype_and_transform(user_curve)


def ask_transtype_and_transform(user_curve):
    print("Available transformations are:")

    # Output available transformations
    for key in user_curve.available_transformations.keys():
        print(f"'{key}'")

    # Includes input validation.
    # Calls input().
    user_input_transformation_key = get_user_input_transformation_key(user_curve)

    # After input validation it should be guaranteed, that user_input_transformation_key is a valid, specific expression for a transformation.

    # Perform the transformation, which will also output the result.
    try:
        user_curve.available_transformations[user_input_transformation_key](print_result=True)

    except Exception as e:
        print_error_and_restart(e, "Transformation")


def print_error_and_restart(e, what_failed):
    print(120 * "-")  # Improves reading comfort.
    print(f"Error: {e}")
    print(f"{what_failed} failed. Please try again with a different curve type or different values (according to error message).")
    print()
    main(print_welcome=False)  # Don't print welcome message.
    sys.exit()


def init_curve(curve_arguments_dic, curves_types, user_input):
    curve_values_tuple = tuple(curve_arguments_dic.values())

    # If var y equals the tuple of the values:
    # user_curve = curves_types[user_input](y[0], y[1], y[2], y[3], y[4])
    # curves_types[user_input] calls the init()/constructor of the selected curve.
    # Initialize the base curve based on user inputs for type and arguments. 
    try:
        user_curve = curves_types[user_input](*curve_values_tuple)

    except Exception as e:
        print_error_and_restart(e, "Init")

    return user_curve


def error_handling_decorator(function):
    # To reuse the error handling of inputs wherever necessary.
    # Avoid copying code.
    def error_handling(*args, **kwargs):
        while True:
            try:

                # Here code that can throw errors.
                result = function(*args, **kwargs)

                return result  # Ends loop.

            except (NameError, SyntaxError):
                # Case: Invalid inputs in eval(str).
                print("Please input a valid integer or mathematical expression in python syntax, e.g. for potency: base**power-subtrahend.")
                continue
            except KeyboardInterrupt:
                print()  # As to not have the new terminal prompt directly after ^C.
                sys.exit()
            # Secure Programming: An attacker should not be able to crash the program.
            except Exception as e:
                print(f"Error: {e}")
                continue

    return error_handling


@error_handling_decorator
def get_user_input_transformation_key(user_curve):
    # Get user input for the key of the transformation.
    user_input_transformation_key = input("Please enter the expression for the curve transformation you want to perform: ")

    # Same input control like first time.
    input_control_key(user_curve.available_transformations, user_input_transformation_key, 6)  # Currently the longest key is 6 chars long: "to ted".

    return user_input_transformation_key


@error_handling_decorator
def eval_parameter_input(key):
    return eval(input(f"Please enter the value for {key}: "))


@error_handling_decorator
def get_user_input_base_key(curves_types):
    # Get user input for the key of the base curve type.
    user_input = input("Please enter the key for the type of curve you want to transform: ")
    input_control_key(curves_types, user_input)
    return user_input


def input_control_key(dic, user_input, maximum_key_length=3):
    # user_input needs input validation.
    # Check if user_input is longer than possible.
    if len(user_input) > maximum_key_length:  # Longest key is maximum_key_length chars long.
        raise ValueError("The input for the key is longer than all keys. Please review your input.")
    # Check if user_input is not a key.
    if user_input not in dic.keys():
        raise ValueError("The input is not a valid key.")


if __name__ == "__main__":
    main(print_welcome=True)  # Print welcome message.

