#include <Std.h>

#include <MpiSupport.h>
#include <Definitions.h>
#include <MathObjects.h>
#include <MathFunctions.h>
#include <IterationCallback.h>

///////////////////////////////////////////////////////////////////////////////

class CExchangeDefinition { // описание одного обмена
public:
    CExchangeDefinition(size_t rank, // ранк того, кому посылаем. откуда == текущий процесс
                        const CMatrixPart &sendPart, // отправляемая часть матрицы
                        const CMatrixPart &recvPart) : // получаемая часть матрицы
            rank(rank),
            sendPart(sendPart),
            recvPart(recvPart) {
        sendBuffer.reserve(sendPart.Size()); // вектор-переменная для обмена
        recvBuffer.reserve(recvPart.Size()); // вектор-переменная для обмена
    }

    const CMatrixPart &SendPart() const { return sendPart; } // getter

    const CMatrixPart &RecvPart() const { return recvPart; } // getter

    void DoExchange(CMatrix &matrix); // асинхронный обмен

    void Wait(CMatrix &matrix); // дождаться обмена

private:
    size_t rank;

    CMatrixPart sendPart; // отправляемая часть матрицы
    MPI_Request sendRequest;
    vector <NumericType> sendBuffer; // вектор-переменная для обмена/ данные запроса на отправку данных в другой процесс

    CMatrixPart recvPart; // получаемая часть матрицы
    MPI_Request recvRequest; // данные запроса на отправку данных в другой процесс
    vector <NumericType> recvBuffer; // вектор-переменная для обмена/ данные запроса на получения данных в другой процесс
};

///////////////////////////////////////////////////////////////////////////////

void CExchangeDefinition::DoExchange(CMatrix &matrix) {
    sendBuffer.clear(); // Очищаем значения, которые посылали в прошлый раз (это вектор)
    for (size_t x = sendPart.BeginX;
         x < sendPart.EndX; x++) { // пушим в буффер для отправки нужную часть матрицы (часть колонки или стобца)
        for (size_t y = sendPart.BeginY; y < sendPart.EndY; y++) {
            sendBuffer.push_back(matrix(x, y));
        }
    }

    MpiCheck( // проверить на 0 и выбросить excepion
            MPI_Isend( // возвращает код ошибки или 0
                    sendBuffer.data(), // адресс начала данных
                    sendBuffer.size(), // кол-во данных
                    MpiNumericType, // тип данных
                    rank, // адресс отправки
                    0, // id сообщения
                    MPI_COMM_WORLD, // коммуникатор
                    &sendRequest) // OUT - "запрос обмена".
            , "MPI_Isend"); // текс exceptionа

    recvBuffer.resize(recvPart.Size());
    MpiCheck( // проверить на 0 и выбросить excepion
            MPI_Irecv( // возвращает код ошибки или 0
                    recvBuffer.data(), // адресс начала данных
                    recvBuffer.size(), // кол-во данных
                    MpiNumericType, // тип данных
                    rank, // адресс отправки
                    0, // id сообщения
                    MPI_COMM_WORLD, // коммуникатор
                    &recvRequest), // OUT - "запрос обмена".
            "MPI_Irecv"); // текс exceptionа
}

void CExchangeDefinition::Wait(CMatrix &matrix) {
    MpiCheck( // проверить на 0 и выбросить excepion
            MPI_Wait( // блокируемся, пока не получим
                    &recvRequest, // переменная, отвечающая за текущий запрос
                     MPI_STATUS_IGNORE), // Mpi_status field будет проигнорирован (передается в Iresv для синхронных операций),
            // иначе можно указатель, куда записывать указать
            "MPI_Wait"); // текс exceptionа

    vector<NumericType>::const_iterator value = recvBuffer.begin(); // Копируем данные в матрицу в текущий процесс
    for (size_t x = recvPart.BeginX; x < recvPart.EndX; x++) {
        for (size_t y = recvPart.BeginY; y < recvPart.EndY; y++) {
            matrix(x, y) = *value;
            ++value;
        }
    }

    MpiCheck( // проверить на 0 и выбросить excepion
            MPI_Wait( // блокируемся, пока не отправим
                    &sendRequest, // переменная, отвечающая за текущий запрос
                    MPI_STATUS_IGNORE),  // Mpi_status field будет проигнорирован (передается в Iresv для синхронных операций),
            // иначе можно указатель, куда записывать указать
            "MPI_Wait"); // текс exceptionа
}

///////////////////////////////////////////////////////////////////////////////

class CExchangeDefinitions : public vector<CExchangeDefinition> { // список обменов
public:
    CExchangeDefinitions() {}

    void Exchange(CMatrix &matrix) { // процедуа выполнения обмена
        for (vector<CExchangeDefinition>::iterator i = begin(); i != end(); ++i) {
            i->DoExchange(matrix); // асинхронный метод обмена
        }
        for (vector<CExchangeDefinition>::iterator i = begin(); i != end(); ++i) {
            i->Wait(matrix); // ждем окончания обмена
        }
    }
};

///////////////////////////////////////////////////////////////////////////////

void GetBeginEndPoints(const size_t numberOfPoints, const size_t numberOfBlocks /* кол-во блоков по абциссе или ординате */,
                       const size_t blockIndex /* текущий номер блока */, size_t &beginPoint,
                       size_t &endPoint) { // Считаем начало и конец отрезка абциссы или ординаты, обрабатываемого процессом
    const size_t objectsPerProcess = numberOfPoints / numberOfBlocks;
    const size_t additionalPoints = numberOfPoints % numberOfBlocks;
    beginPoint = objectsPerProcess * blockIndex + min(blockIndex, additionalPoints);
    endPoint = beginPoint + objectsPerProcess;
    if (blockIndex < additionalPoints) {
        endPoint++;
    }
}

///////////////////////////////////////////////////////////////////////////////

NumericType TotalError(const CMatrix &p, const CUniformGrid &grid) { // Считаем невязку
    NumericType squares = 0;
    for (size_t x = 1; x < p.SizeX() - 1; x++) {
        for (size_t y = 1; y < p.SizeY() - 1; y++) { // Считаем сумму квадратов
            squares += (Phi(grid.X[x], grid.Y[y]) - p(x, y)) * (Phi(grid.X[x], grid.Y[y]) - p(x, y));
        }
    }
    return static_cast<NumericType> (pow(squares, 0.5));
}

///////////////////////////////////////////////////////////////////////////////

// Вывод результатов (матрицы) в файл
void DumpMatrix(const CMatrix &matrix, const CUniformGrid &grid, ostream &output) {
    for (size_t x = 0; x < matrix.SizeX(); x++) {
        for (size_t y = 0; y < matrix.SizeY(); y++) {
            output << grid.X[x] << '\t' << grid.Y[y] << '\t' << matrix(x, y) << endl;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

class CProgram {
public:
    static void Run(size_t pointsX, size_t pointsY, const CArea &area,
                    IIterationCallback &callback, const string &dumpFilename);

private:
    const size_t numberOfProcesses;
    const size_t rank;
    const size_t pointsX; // число узлов сетки
    const size_t pointsY;
    size_t processesX; // число MPI процессов "обрабатывающих оси"
    size_t processesY;
    size_t rankX; // Порядковый номер прямоугольника
    size_t rankY;
    size_t beginX; // Границы обрабатываемого прямоугольника
    size_t endX;
    size_t beginY;
    size_t endY;
    CExchangeDefinitions exchangeDefinitions; // С кем и чем обменивается процесс
    CUniformGrid grid;
    CMatrix p; // Приближение
    CMatrix r; // Направление движения к следующему приближжению на 1 итерации
    CMatrix g; // Направление движения к следующему приближжению
    NumericType difference; // Невязка
    NumericType difference_2; // Сумма квадратов разниц (для AllReduceDifference)

    CProgram(size_t pointsX, size_t pointsY, const CArea &area);

    bool hasLeftNeighbor() const { return (rankX > 0); }

    bool hasRightNeighbor() const { return (rankX < (processesX - 1)); }

    bool hasTopNeighbor() const { return (rankY > 0); }

    bool hasBottomNeighbor() const { return (rankY < (processesY - 1)); }

    size_t rankByXY(size_t x, size_t y) const { return (y * processesX + x); }

    void setProcessXY(); // узнаем, сколько процессов будет по абциссе, сколько по ординате, деля их число на 2...

    void setExchangeDefinitions(); // Заполняем список соседей с которыми будем обмениваться данными и описываем сами данные.

    void allReduceFraction(CFraction &fraction);

    void allReduceDifference();

    void iteration0(); // итерация 0 == инициализация матрицы

    void iteration1(); // итерация 1, выполняется по отдельной формуле

    void iteration2(); // остальные итерации, для ускорения, см. методичку
};

///////////////////////////////////////////////////////////////////////////////

void CProgram::Run(size_t pointsX, size_t pointsY, const CArea &area,
                   IIterationCallback &callback, const string &dumpFilename = "") {
    CProgram program(pointsX, pointsY, area); // Конструктор запускаем

    // Выполняем нулевую итерацию (инициализацию).
    if (!callback.BeginIteration()) {
        return;
    }
    program.iteration0(); // Заполняем границы, если границы общей области принадлежат области, обрабатываемой процессом
    callback.EndIteration(program.difference); // Значение difference по умолчанию задается в конструкторе program

    // Выполняем первую итерацию.
    if (!callback.BeginIteration()) {
        return;
    }
    program.iteration1();
    callback.EndIteration(program.difference);

    // Выполняем остальные итерации.
    while (callback.BeginIteration()) { // проверяем невязку
        program.iteration2(); // выполняем итерацию
        callback.EndIteration(program.difference); // проставляем невязку и логгируем итерацию
    }

    char num[5];
    snprintf(num, 5, "%d", (int) CMpiSupport::Rank()); // данные текущего процесса записываются в файл с именем +  mpi-ранк процесса
    ofstream outputFile((dumpFilename + string(num)).c_str());
    DumpMatrix(program.p, program.grid, outputFile); // выводим нашу матрицу
}

CProgram::CProgram(size_t pointsX, size_t pointsY, const CArea &area) :
        numberOfProcesses(CMpiSupport::NumberOfProccess()),
        rank(CMpiSupport::Rank()),
        pointsX(pointsX), pointsY(pointsY),
        difference(numeric_limits<NumericType>::max()) {
    setProcessXY(); // сколько процессов "по горизонтали и вертикали" (для 8 4:2)
    rankX = rank % processesX; // какую часть обрабатывает этот процесс
    rankY = rank / processesX;
    GetBeginEndPoints(pointsX, processesX, rankX, beginX,
                      endX); // Считаем начало и конец отрезка, обрабатываемого процессом
    GetBeginEndPoints(pointsY, processesY, rankY, beginY,
                      endY); // Считаем начало и конец отрезка, обрабатываемого процессом

    if (hasLeftNeighbor()) { // Корректируем концы, чтобы было с "заездом" на чужую территорию
        beginX--;
    }
    if (hasRightNeighbor()) {
        endX++;
    }
    if (hasTopNeighbor()) {
        beginY--;
    }
    if (hasBottomNeighbor()) {
        endY++;
    }

    // Инициализируем grid.
    grid.X.PartInit(area.X0, area.Xn, pointsX, beginX, endX); // У каждого процесса свой грид
    grid.Y.PartInit(area.Y0, area.Yn, pointsY, beginY, endY);

    // Заполняем список соседей с которыми будем обмениваться данными и описываем сами данные.
    setExchangeDefinitions();
}

void CProgram::setProcessXY() { // узнаем, сколько процессов будет по абциссе, сколько по ординате, деля их число на 2...
    size_t power = 0;
    {
        size_t i = 1;
        while (i < numberOfProcesses) {
            i *= 2;
            power++;
        }
        if (i != numberOfProcesses) {
            throw CException("The number of processes must be power of 2.");
        }
    }

    float pX = static_cast<float>( pointsX );
    float pY = static_cast<float>( pointsY );

    size_t powerX = 0;
    size_t powerY = 0;
    for (size_t i = 0; i < power; i++) {
        if (pX > pY) {
            pX = pX / 2;
            powerX++;
        } else {
            pY = pY / 2;
            powerY++;
        }
    }

    processesX = 1 << powerX;
    processesY = 1 << powerY;
}

void CProgram::setExchangeDefinitions() {
    if (hasLeftNeighbor()) { // проверка, не крайний ли наш блок
        exchangeDefinitions.push_back(CExchangeDefinition(
                rankByXY(rankX - 1, rankY), // Устанавливаем ранк соседа, с которым будем обмениваться
                grid.Column(1, 1 /* decreaseTop */, 1 /* decreaseBottom */ ),
                grid.Column(0, 1 /* decreaseTop */, 1 /* decreaseBottom */ )));
    }
    if (hasRightNeighbor()) { // проверка, не крайний ли наш блок
        exchangeDefinitions.push_back(CExchangeDefinition(
                rankByXY(rankX + 1, rankY), // Устанавливаем ранк соседа, с которым будем обмениваться
                grid.Column(grid.X.Size() - 2, 1 /* decreaseTop */, 1 /* decreaseBottom */ ),
                grid.Column(grid.X.Size() - 1, 1 /* decreaseTop */, 1 /* decreaseBottom */ )));
    }
    if (hasTopNeighbor()) { // проверка, не крайний ли наш блок
        exchangeDefinitions.push_back(CExchangeDefinition(
                rankByXY(rankX, rankY - 1), // Устанавливаем ранк соседа, с которым будем обмениваться
                grid.Row(1, 1 /* decreaseLeft */, 1 /* decreaseRight */ ),
                grid.Row(0, 1 /* decreaseLeft */, 1 /* decreaseRight */ )));
    }
    if (hasBottomNeighbor()) { // проверка, не крайний ли наш блок
        exchangeDefinitions.push_back(CExchangeDefinition(
                rankByXY(rankX, rankY + 1), // Устанавливаем ранк соседа, с которым будем обмениваться
                grid.Row(grid.Y.Size() - 2, 1 /* decreaseLeft */, 1 /* decreaseRight */ ),
                grid.Row(grid.Y.Size() - 1, 1 /* decreaseLeft */, 1 /* decreaseRight */ )));
    }
}

void CProgram::allReduceFraction(CFraction &fraction) {
    NumericType buffer[2] = {fraction.Numerator, fraction.Denominator};
    MpiCheck( // проверяем на MPI_SUCCESS == 0
            MPI_Allreduce(MPI_IN_PLACE, // input buffer == output buffer
                          buffer, // данные
                          2, // размер
                          MpiNumericType, // тип
                          MPI_SUM, // операция
                          MPI_COMM_WORLD), // коммуникатор
            "MPI_Allreduce" // текст ошибки
    );
    fraction.Numerator = buffer[0]; // числитель
    fraction.Denominator = buffer[1]; // знаменатель
}

void CProgram::allReduceDifference() {
    NumericType buffer = difference_2;
    MpiCheck( // проверяем на MPI_SUCCESS == 0
            MPI_Allreduce(MPI_IN_PLACE, // input buffer == output buffer
                          &buffer, // адресс переменной, с данными запроса-ответа
                          1, // размер
                          MpiNumericType, // тип
                          MPI_SUM, // суммируем квадраты
                          MPI_COMM_WORLD),
            "MPI_Allreduce" // текст ошибки
    );
    difference = static_cast<NumericType> (pow(buffer, 0.5)); // считаем общую невязку
}

void CProgram::iteration0() {
    // Заполняем границы, если границы общей области принадлежат области, обрабатываемой процессом
    p.Init(grid.X.Size(), grid.Y.Size());

    if (!hasLeftNeighbor()) {
        for (size_t y = 0; y < p.SizeY(); y++) {
            p(0, y) = Phi(grid.X[0], grid.Y[y]); // Phi - граничная функция
        }
    }
    if (!hasRightNeighbor()) {
        const size_t left = p.SizeX() - 1;
        for (size_t y = 0; y < p.SizeY(); y++) {
            p(left, y) = Phi(grid.X[left], grid.Y[y]); // Phi - граничная функция
        }
    }
    if (!hasTopNeighbor()) {
        for (size_t x = 0; x < p.SizeX(); x++) {
            p(x, 0) = Phi(grid.X[x], grid.Y[0]); // Phi - граничная функция
        }
    }
    if (!hasBottomNeighbor()) {
        const size_t bottom = p.SizeY() - 1;
        for (size_t x = 0; x < p.SizeX(); x++) {
            p(x, bottom) = Phi(grid.X[x], grid.Y[bottom]); // Phi - граничная функция
        }
    }
}

void CProgram::iteration1() {
    r.Init(grid.X.Size(), grid.Y.Size());

    CalcR(p, grid, r);
    exchangeDefinitions.Exchange(r);

    CFraction tau = CalcTau(r, r, grid);
    allReduceFraction(tau);

    difference = CalcP(r, tau.Value(), p);
    difference_2 = CalcP_2(r, tau.Value(), p);
    allReduceDifference();

    g = r;
}

void CProgram::iteration2() {
    exchangeDefinitions.Exchange(p);

    CalcR(p, grid, r);
    exchangeDefinitions.Exchange(r);

    CFraction alpha = CalcAlpha(r, g, grid);
    allReduceFraction(alpha);

    CalcG(r, alpha.Value(), g);
    exchangeDefinitions.Exchange(g);

    CFraction tau = CalcTau(r, g, grid);
    allReduceFraction(tau);

    difference = CalcP(g, tau.Value(), p);
    difference_2 = CalcP_2(g, tau.Value(), p);
    allReduceDifference();
}

///////////////////////////////////////////////////////////////////////////////

// Последовательная реализация.
void Serial(const size_t pointsX, const size_t pointsY, const CArea &area,
            IIterationCallback &callback, const string &dumpFilename = "") {
    // Инициализируем grid.
    CUniformGrid grid;
    grid.X.Init(area.X0, area.Xn, pointsX); // MathObjects.cpp -> PartInit(0 3 100 0 100)
    grid.Y.Init(area.Y0, area.Yn, pointsY);

    CMatrix p(grid.X.Size(), grid.Y.Size()); // create empty matrixes
    CMatrix r(grid.X.Size(), grid.Y.Size());

    NumericType difference = numeric_limits<NumericType>::max(); // max NumericType

    // Выполняем нулевую итерацию (инициализацию).
    if (!callback.BeginIteration()) { // if diff-eps and max number of iteration bad
        return;
    }

    for (size_t x = 0; x < p.SizeX(); x++) {
        p(x, 0) = Phi(grid.X[x], grid.Y[0]); // border values
        p(x, p.SizeY() - 1) = Phi(grid.X[x], grid.Y[p.SizeY() - 1]); // border values
    }
    for (size_t y = 1; y < p.SizeY() - 1; y++) {
        p(0, y) = Phi(grid.X[0], grid.Y[y]);  // border values
        p(p.SizeX() - 1, y) = Phi(grid.X[p.SizeX() - 1], grid.Y[y]);  // border values
    }
    callback.EndIteration(difference); // with max NumericType

    // Выполняем первую итерацию.
    if (!callback.BeginIteration()) {
        return;
    }
    {
        CalcR(p, grid, r); // Cчитаем невязку r в неграничных точках
        const CFraction tau = CalcTau(r, r, grid); // считаем tau_1
        difference = CalcP(r, tau.Value(), p); // Вычисление значений pij во внутренних точках, возвращается норма.
    }
    callback.EndIteration(difference);

    CMatrix g(r);
    // Выполняем остальные итерации.
    while (callback.BeginIteration()) { // выйдем из цикла, когда достигнем eps
        CalcR(p, grid, r); // Cчитаем невязку r в неграничных точках
        const CFraction alpha = CalcAlpha(r, g, grid); // параметр скорейшего спуска
        CalcG(r, alpha.Value(), g); // считаем направление
        const CFraction tau = CalcTau(r, g, grid); // считаем tau_k
        difference = CalcP(g, tau.Value(), p); // Вычисление значений pij во внутренних точках, возвращается норма.

        callback.EndIteration(difference);
    }

    if (!dumpFilename.empty()) { // 3 аргумент - вывод результата
        cout << "Total error: " << TotalError(p, grid) << endl;
        ofstream outputFile(dumpFilename.c_str());
        DumpMatrix(p, grid, outputFile);
    }
}

///////////////////////////////////////////////////////////////////////////////

void ParseArguments(const int argc, const char *const argv[],
                    size_t &pointsX, size_t &pointsY, string &dumpFilename) { // read arguments
    if (argc < 3 || argc > 4) {
        throw CException("too few arguments\n"
                                 "Usage: dirch POINTS_X POINTS_Y [DUMP_FILENAME]");
    }

    pointsX = strtoul(argv[1], 0, 10);
    pointsY = strtoul(argv[2], 0, 10);

    if (pointsX == 0 || pointsY == 0) {
        throw CException("invalid format of arguments\n"
                                 "Usage: dirch POINTS_X POINTS_Y [DUMP_FILENAME]");
    }

    if (argc == 4) {
        dumpFilename = argv[3];
    }
}

void Main(const int argc, const char *const argv[]) {
    double programTime = 0.0; // default timer time
    {
        CMpiTimer timer(programTime); // initialize timer

        // Используем только потоки ввода вывода iostream,
        // поэтому отключаем синхронизацию ввода вывода со стандартной библиотекой C.
        ios::sync_with_stdio(false);

        size_t pointsX;
        size_t pointsY;
        string dumpFilename;
        ParseArguments(argc, argv, pointsX, pointsY, dumpFilename); // read arguments

        auto_ptr <IIterationCallback> callback(new CSimpleIterationCallback);
        if (CMpiSupport::Rank() == 0) { // if main mpi process
            callback.reset(new CIterationCallback(cout, 0)); // destruct and create new
        }

        if (CMpiSupport::NumberOfProccess() == 1) { // only one process
            Serial(pointsX, pointsY, Area, *callback, dumpFilename);
        } else { // more then one process
            CProgram::Run(pointsX, pointsY, Area, *callback, dumpFilename);
        }
    }
    cout << "(" << CMpiSupport::Rank() << ") Time: " << programTime << endl;
}

// Точка старта
int main(int argc, char **argv) // start point
{
    try {
        CMpiSupport::Initialize(&argc, &argv); // initialize mpi
        Main(argc, argv); // main function
        CMpiSupport::Finalize(); // finalize mpi
    } catch (exception &e) {
        cerr << "Error: " << e.what() << endl;
        CMpiSupport::Abort(1); // abort or throw
        return 1;
    } catch (...) {
        cerr << "Unknown error!" << endl;
        CMpiSupport::Abort(2); // abort or throw
        return 2;
    }

    return 0;
}
