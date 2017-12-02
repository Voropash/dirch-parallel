#include <Std.h>
#include <Errors.h>
#include <MpiSupport.h>

///////////////////////////////////////////////////////////////////////////////

string CMpiException::makeErrorText() const { // Вызывается в конструкторе exception
    ostringstream oss;
    oss << "MPI function '" << FunctionName() << "'"
        << " has failed with code '" << ErrorCode() << "'.";
    return oss.str();
}

///////////////////////////////////////////////////////////////////////////////

void MpiCheck(const int mpiResult, const string &mpiFunctionName) { // check for 0 and throw...
    if (mpiResult != MPI_SUCCESS) {  // MPI_SUCCESS == 0
        throw CMpiException(mpiResult, mpiFunctionName);
    }
}

///////////////////////////////////////////////////////////////////////////////

bool CMpiSupport::initialized = false; // default value
size_t CMpiSupport::rank = 0;
size_t CMpiSupport::numberOfProccess = 0;

void CMpiSupport::Initialize(int *argc, char ***argv) {
    if (Initialized()) { // it should be first initialization
        throw CException("MPI was already initialized!");
    }
    MpiCheck(MPI_Init(argc, argv), "MPI_Init"); // init and throw if error
    int tmp;
    MpiCheck(MPI_Comm_rank(MPI_COMM_WORLD, &tmp), "MPI_Comm_rank"); // get rank and throw if error
    rank = static_cast<size_t>( tmp );
    MpiCheck(MPI_Comm_size(MPI_COMM_WORLD, &tmp), "MPI_Comm_size"); // get number of process and throw if error
    numberOfProccess = static_cast<size_t>( tmp );
    initialized = true;
}

void CMpiSupport::Finalize() {
    checkInitialized(); // throw if not initialized
    MPI_Finalize();
}

void CMpiSupport::Abort(int code) { // abort processing
    if (Initialized()) {
        MPI_Abort(MPI_COMM_WORLD, code); // send "please abort"
    }
}

void CMpiSupport::checkInitialized() { // check and throw
    if (!Initialized()) {
        throw CException("MPI was not initialized yet!");
    }
}

void CMpiSupport::Barrier() // make the barrier, throw exc in case of error
{
    MpiCheck(MPI_Barrier(MPI_COMM_WORLD /* кол-во процессов */),
             "MPI_Barrier" /* function name for loggin in stacktrace*/ );
}

///////////////////////////////////////////////////////////////////////////////
