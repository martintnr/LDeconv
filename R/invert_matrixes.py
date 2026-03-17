import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from numbers import Real

import numpy as np
from scipy.sparse import csc_matrix, save_npz
from scipy.sparse.linalg import eigsh
from sklearn.utils.extmath import randomized_svd


# Perform truncated SVD and compute the pseudo-inverse of the matrix
def fast_truncated_svd_pinv(submatrix, k, random_state=42):
    U, S, Vt = randomized_svd(submatrix, n_components=k, random_state=random_state)
    S_inv = np.diag(1.0 / S)
    return Vt.T @ S_inv @ U.T


# Choose the number of components based on the amount of variance explained
def choose_k_with_eigenvalues(submatrix, variance_threshold=0.99):
    eigenvalues, _ = eigsh(submatrix, k=submatrix.shape[0] - 1, which="LM")
    variance_explained_ratio = np.cumsum(np.sort(eigenvalues)[::-1]) / np.sum(eigenvalues)
    k = np.argmax(variance_explained_ratio >= variance_threshold) + 1
    return eigenvalues, k


def build_refined_triplets(categories):
    base_blocks = {}

    for i, cat in enumerate(categories):
        if isinstance(cat, Real) and not isinstance(cat, bool):
            if float(cat).is_integer():
                cat = int(cat)
        elif isinstance(cat, str):
            try:
                x = float(cat)
                if x.is_integer():
                    cat = int(x)
            except ValueError:
                pass

        base_blocks.setdefault(cat, []).append(i)

    sorted_keys = sorted(base_blocks.keys())
    refined_blocks = {}

    for i, key in enumerate(sorted_keys):
        merged_indices = []

        if i > 0:
            merged_indices.extend(base_blocks[sorted_keys[i - 1]])

        merged_indices.extend(base_blocks[key])

        if i < len(sorted_keys) - 1:
            merged_indices.extend(base_blocks[sorted_keys[i + 1]])

        refined_blocks[key] = merged_indices

    return refined_blocks


# Invert a single LD block and save it as a sparse matrix file
def invert_and_save_block(args):
    block_id, block, variance_threshold, zero_tol, output_dir, chromosome = args

    try:
        out_path = os.path.join(output_dir, f"chr{chromosome}_LD_inv_BD_{block_id}.npz")
        if os.path.exists(out_path):
            return block_id

        if block.shape == (1, 1):
            a = block[0, 0]
            a = a if isinstance(a, (int, float)) else float(a)
            inv = [[0.0]] if abs(a) <= zero_tol else [[1.0 / a]]
            save_npz(out_path, csc_matrix(inv))
            return block_id

        eigenvalues, k = choose_k_with_eigenvalues(block, variance_threshold)
        inverted_block = fast_truncated_svd_pinv(block, k)
        save_npz(out_path, csc_matrix(inverted_block))
        return block_id

    except Exception as e:
        print(f"Error inverting block {block_id} (chr{chromosome}): {e}")
        return None


# Run block inversions
def invert_ld_blocks_parallel_and_save(
    LD_blocks_diag,
    variance_threshold,
    zero_tol,
    output_dir,
    chromosome,
):
    results = []

    for bid in sorted(LD_blocks_diag.keys()):
        block = LD_blocks_diag[bid]
        res = invert_and_save_block(
            (bid, block, variance_threshold, zero_tol, output_dir, chromosome)
        )
        if res is not None:
            results.append(res)

    return results


# ---- Main execution starts here ----

input_dir = "/nas/pleiomap/Data/UKBB_LD_sparse_mat_R"
output_dir = "/nas/pleiomap/Data/UKBB_inverted_LD_AdjClust_099_10KAbsR_SW"
os.makedirs(output_dir, exist_ok=True)


def process_chromosomes(
    input_dir,
    output_dir,
    variance_threshold=0.99,
    zero_tol=1e-12,
    start=1,
    stop=23,
    step=1,
    parallel=True,
    max_workers=4,
):
    def _process_one(chromosome: int):
        print(f"\nProcessing chromosome {chromosome}...")

        file_name = f"chr{chromosome}_sparse_matrix.npz"
        ld_path = os.path.join(input_dir, file_name)

        ld = np.load(ld_path)
        LD = csc_matrix((ld["data"], ld["indices"], ld["indptr"]), shape=ld["shape"])

        rlf = f"chr{chromosome}_refined_block_labels.npy"
        rlf_path = os.path.join(output_dir, rlf)

        refined_labels = np.load(rlf_path)
        refined_triplets = build_refined_triplets(refined_labels)

        LD_blocks_diag = {
            bid: LD[np.ix_(indices, indices)]
            for bid, indices in refined_triplets.items()
        }

        invert_ld_blocks_parallel_and_save(
            LD_blocks_diag,
            variance_threshold=variance_threshold,
            zero_tol=zero_tol,
            output_dir=output_dir,
            chromosome=chromosome,
        )

        return chromosome

    chromosomes = list(range(start, stop, step))

    if parallel:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(_process_one, chrom): chrom for chrom in chromosomes}
            for fut in as_completed(futures):
                chrom = futures[fut]
                try:
                    fut.result()
                except Exception as e:
                    print(f"Chromosome {chrom} failed: {e}")
    else:
        for chrom in chromosomes:
            try:
                _process_one(chrom)
            except Exception as e:
                print(f"Chromosome {chrom} failed: {e}")


process_chromosomes(
    input_dir=input_dir,
    output_dir=output_dir,
    variance_threshold=0.99,
    zero_tol=1e-2,
    parallel=True,
    max_workers=2,
)
