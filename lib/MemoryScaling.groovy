// Shared task.attempt -> memory (GB) retry series for OOM-prone processes across
// the pipeline -- see params.max_memory in nextflow.config. Three series depending
// on what a process realistically needs on a healthy first attempt: each follows a
// short hand-tuned ramp, then doubles every attempt after that, but never asks for
// more than params.max_memory however many attempts a process' own maxRetries
// allows.
//
//   SERIES_4  (processes that need ~4GB on attempt 1, e.g. pileup_variants):
//             4, 8, 12, 16, 24, 32, 64, 128, 256, ...
//   SERIES_8  (processes that need ~8GB on attempt 1, e.g. refine_with_sv):
//             8, 16, 24, 32, 64, 128, 256, ...
//   SERIES_16 (processes that need ~16GB on attempt 1, e.g. run_tapes):
//             16, 24, 32, 64, 128, 256, ...
class MemoryScaling {

    static final List<Integer> SERIES_4  = [4, 8, 12, 16, 24, 32]
    static final List<Integer> SERIES_8  = [8, 16, 24, 32]
    static final List<Integer> SERIES_16 = [16, 24, 32]

    // memory to request on this task.attempt (1-based), as a Nextflow memory
    // directive string (e.g. "16 GB") -- never exceeds maxMemoryGB.
    static String forAttempt(List<Integer> baseGB, int attempt, int maxMemoryGB) {
        long gb = gbForAttempt(baseGB, attempt)
        if (gb > maxMemoryGB) {
            gb = maxMemoryGB
        }
        return "${gb} GB"
    }

    // how many maxRetries (i.e. attempts - 1) are needed for forAttempt() to reach
    // (or first meet/exceed) maxMemoryGB -- lets each process' maxRetries directive
    // stay correct if params.max_memory is overridden away from its 128GB default.
    static int retriesNeeded(List<Integer> baseGB, int maxMemoryGB) {
        int attempt = 1
        while (gbForAttempt(baseGB, attempt) < maxMemoryGB) {
            attempt++
        }
        return attempt - 1
    }

    private static long gbForAttempt(List<Integer> baseGB, int attempt) {
        int idx = attempt - 1
        if (idx < baseGB.size()) {
            return baseGB[idx]
        }
        int doublings = idx - baseGB.size() + 1
        return baseGB[-1] * (1L << doublings)
    }
}
