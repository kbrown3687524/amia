import os
import argparse
import requests
import logging
from datetime import datetime
from multiprocessing import Process


class AutoAllosteric:

    def submit_to_passer(pdb_path, model):
        url = "https://passer.smu.edu/api"
        with open(pdb_path, 'rb') as pdb_file:
            files = {'pdbFile': pdb_file}
            data = {"model": model, "format": "json"}
            try:
                logging.info(f"📤 Submitting {pdb_path} to PASSer with model '{model}'")
                response = requests.post(url, files=files, data=data)
                response.raise_for_status()
                logging.info(f"✅ Received response for {os.path.basename(pdb_path)} ({model})")
                return response.json()
            except Exception as e:
                logging.error(f"❌ Failed submission for {pdb_path} ({model}): {e}")
                return {"error": str(e)}

    @staticmethod
    def save_passer_text_results(pdb_dir, output_text):
        models = ["ensemble", "automl", "rank"]
        with open(output_text, "w") as out:
            out.write("PASSer batch submission results\n")
            out.write("=" * 40 + "\n\n")

        for filename in os.listdir(pdb_dir):
            if filename.endswith("_auto.pdb"):
                pdb_path = os.path.join(pdb_dir, filename)
                for model in models:
                    result = AutoAllosteric.submit_to_passer(pdb_path, model)
                    with open(output_text, "a") as out:
                        out.write(f"\n=== Results for {filename} ({model}) ===\n")
                        out.write(str(result) + "\n")
                    logging.info(f"📝 Written results for {filename} ({model}) to {output_text}")

    @staticmethod
    def download_passer_zip_results(pdb_dir, output_dir):
        os.makedirs(output_dir, exist_ok=True)
        models = ["ensemble", "automl", "rank"]

        for filename in os.listdir(pdb_dir):
            if filename.endswith("_auto.pdb"):
                pdb_path = os.path.join(pdb_dir, filename)
                base_name = os.path.splitext(filename)[0]
                for model in models:
                    try:
                        url = "https://passer.smu.edu/api"
                        with open(pdb_path, 'rb') as pdb_file:
                            files = {'pdbFile': pdb_file}
                            data = {"model": model, "format": "zip"}
                            logging.info(f"📦 Requesting ZIP for {filename} ({model})")
                            response = requests.post(url, files=files, data=data)
                            response.raise_for_status()
                            zip_name = f"{base_name}_{model}.zip"
                            zip_path = os.path.join(output_dir, zip_name)
                            with open(zip_path, 'wb') as f:
                                f.write(response.content)
                            logging.info(f"✅ Saved ZIP to {zip_path}")
                    except Exception as e:
                        logging.error(f"❌ Failed ZIP for {filename} ({model}): {e}")

    @staticmethod
    def process_individual_pdb(pdb_file, output_dir, output_text):
        models = ["ensemble", "automl", "rank"]
        os.makedirs(output_dir, exist_ok=True)
        with open(output_text, "a") as out:
            out.write(f"Results for {os.path.basename(pdb_file)}\n")
            out.write("=" * 40 + "\n\n")

        for model in models:
            result = AutoAllosteric.submit_to_passer(pdb_file, model)
            with open(output_text, "a") as out:
                out.write(f"\n=== Model: {model} ===\n")
                out.write(str(result) + "\n")
            logging.info(f"📝 Added {model} results to {output_text}")

            try:
                with open(pdb_file, 'rb') as pdb_f:
                    files = {'pdbFile': pdb_f}
                    data = {"model": model, "format": "zip"}
                    response = requests.post("https://passer.smu.edu/api", files=files, data=data)
                    response.raise_for_status()
                    zip_name = f"{os.path.splitext(os.path.basename(pdb_file))[0]}_{model}.zip"
                    zip_path = os.path.join(output_dir, zip_name)
                    with open(zip_path, 'wb') as f:
                        f.write(response.content)
                    logging.info(f"✅ Saved ZIP to {zip_path}")
            except Exception as e:
                logging.error(f"❌ Failed to download ZIP for {model}: {e}")

    @staticmethod
    def process_output_file(text_output, html_output):
        try:
            with open(text_output, 'r') as f:
                lines = f.readlines()

            with open(html_output, 'w') as html:
                html.write("<html><body><h1>PASSer Summary</h1><pre>\n")
                html.writelines(lines)
                html.write("</pre></body></html>\n")

            logging.info(f"🖥️ Generated HTML summary: {html_output}")
        except Exception as e:
            logging.error(f"❌ Failed to generate HTML summary: {e}")


def setup_logger(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, f"passer_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")

    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
    logging.getLogger().addHandler(console)

    logging.info("🚀 Logger initialized.")
    logging.info(f"📁 Logs will be saved to {log_file}")


def main():
    parser = argparse.ArgumentParser(description="Submit PDB files to PASSer and process results.")
    parser.add_argument('--pdb_file', required=False, help='Optional path to a single PDB file to process.')
    parser.add_argument('--pdb_dir', required=False, help='Directory with input PDB files (ending in _auto.pdb).')
    parser.add_argument('--output_dir', required=True, help='Directory to save PASSer ZIP output files.')
    parser.add_argument('--text_output', required=True, help='Text file to save PASSer JSON/text results.')
    parser.add_argument('--html_output', required=True, help='HTML file for pocket score summary.')
    args = parser.parse_args()

    setup_logger(args.output_dir)

    if not args.pdb_file and not args.pdb_dir:
        parser.error("Either --pdb_dir or --pdb_file must be provided.")

    logging.info("🔄 PASSer submission workflow started.")

    if args.pdb_dir:
        logging.info(f"📂 Processing directory: {args.pdb_dir}")
        p1 = Process(target=AutoAllosteric.save_passer_text_results, args=(args.pdb_dir, args.text_output))
        p2 = Process(target=AutoAllosteric.download_passer_zip_results, args=(args.pdb_dir, args.output_dir))
        p1.start()
        p2.start()
        p1.join()
        p2.join()

    if args.pdb_file:
        logging.info(f"📄 Processing single PDB file: {args.pdb_file}")
        AutoAllosteric.process_individual_pdb(args.pdb_file, args.output_dir, args.text_output)

    AutoAllosteric.process_output_file(args.text_output, args.html_output)

    logging.info("✅ Workflow complete.")


if __name__ == "__main__":
    main()