import argparse
from gem.pipeline import run_all

def main():
    parser = argparse.ArgumentParser(
        description="GEM: Genetic Exchange Model CLI"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Define `run-all` subcommand with default values shown in help
    run_parser = subparsers.add_parser(
        "run-all",
        help="Run the full GEM pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    run_parser.add_argument("--target", required=True, help="FASTA file of target gene sequences")
    run_parser.add_argument("--known", required=True, help="FASTA file of known host sequences")
    run_parser.add_argument("--novel", required=True, help="FASTA file of novel host sequences to be linked")
    run_parser.add_argument("--email", required=True, help="Email for NCBI Entrez queries")
    run_parser.add_argument("--min-len", type=int, default=5000, help="Minimum length to keep known sequences")
    run_parser.add_argument("--segment-size", type=int, default=5000, help="Up/downstream window size in bp")
    run_parser.add_argument("--d-range", nargs=3, type=int, default=[0, 12000, 2000],
                            help="Three space-separated integers: start, end, and step of expansion distances")
    run_parser.add_argument("--coverage-threshold", type=int, default=4000, help="Minimum total BLAST alignment length")
    run_parser.add_argument("--identity-threshold", type=float, default=80.0, help="Minimum BLAST identity percentage")
    run_parser.add_argument("--evalue-threshold", type=float, default=1e-3, help="Maximum acceptable BLAST e-value")
    run_parser.add_argument("--threads", type=int, default=1, help="Number of threads for BLAST")

    args = parser.parse_args()

    if args.command == "run-all":
        run_all(
            target=args.target,
            known=args.known,
            novel=args.novel,
            email=args.email,
            min_len=args.min_len,
            segment_size=args.segment_size,
            d_range=tuple(args.d_range),
            coverage_threshold=args.coverage_threshold,
            identity_threshold=args.identity_threshold,
            evalue_threshold=args.evalue_threshold,
            threads=args.threads
        )

if __name__ == "__main__":
    main()
