#!/usr/bin/env python3
import argparse

def fill_template(template_path, output_path, replacements):
    """Replace placeholder keywords in a template XML file."""
    with open(template_path, "r") as f:
        text = f.read()

    # Perform all replacements (exact string match)
    for key, value in replacements.items():
        text = text.replace(key, str(value))

    with open(output_path, "w") as f:
        f.write(text)

    print(f"Written filled XML to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Fill i-PI XML template.")
    parser.add_argument(
        "--template", 
        required=True, 
        help="Path to XML template"
        )
    parser.add_argument(
        "--output", 
        required=True, 
        help="Output XML file"
        )

    # Add arguments for all your CAPS placeholders
    parser.add_argument(
        "--filename", 
        required=True,
        help="Output filename prefix"
        )
    parser.add_argument(
        "--model_path", 
        required=True,
        help="Path to the ML model"
    )
    parser.add_argument(
        "--inputfile", 
        required=True,
        help="Input structure file"
    )
    parser.add_argument(
        "--nsteps", 
        default=20000,
        help="Number of MD steps"
        )
    parser.add_argument(
        "--seed",
        default=12345,
        help="Random seed for PRNG"
    )
    parser.add_argument(
        "--temperature", 
        default=300,
        help="Simulation temperature in K"
    )
    parser.add_argument(
        "--pressure", 
        default=1.0,
        help="Simulation pressure in GPa"
    )

    args = parser.parse_args()

    # Mapping from placeholder → provided value
    repl = {
        "FILENAME": args.filename,
        "NSTEPS": args.nsteps,
        "SEED": args.seed,
        "MODEL_PATH": args.model_path,
        "INPUTFILE": args.inputfile,
        "TEMPERATURE": args.temperature,
        "PRESSURE": args.pressure,
    }

    fill_template(args.template, args.output, repl)


if __name__ == "__main__":
    main()
