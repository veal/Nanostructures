#ifndef TASK_H
#define TASK_H

class Task
{
public:
    Task();
    virtual void operator()() = 0;
};

#endif // TASK_H
