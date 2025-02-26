
#include "uri.h"

static char weekday_name[][20] = {"Sunday",   "Monday", "Tuesday", "Wednesday",
                                  "Thursday", "Friday", "Saturday"};

static char ab_weekday_name[][10] = {"Sun", "Mon", "Tue", "Wed",
                                     "Thu", "Fri", "Sat"};

static char month_name[][20] = {
    "January", "February", "March",     "April",   "May",      "June",
    "July",    "August",   "September", "October", "November", "December"};

static char ab_month_name[][10] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

static const unsigned short int __mon_yday[2][13] = {
    /* Normal years.  */
    {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
    /* Leap years.  */
    {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}};

#define IS_LEAP(year) \
    (0 == (year) % 4 && (0 != (year) % 100 || 0 == (year) % 400))

/* Compute the day of the week.  */
int dayOfWeek(uritm *tm)
{
    /* January 1st 1970 was a Thursday (= 4). */

    int corr_year = tm->year - (tm->month < 2);
    int wday =
        (-473 + (365 * (tm->year - 1970)) + (corr_year / 4) -
         ((corr_year / 4) / 25) + ((corr_year / 4) % 25 < 0) +
         (((corr_year / 4) / 25) / 4) + __mon_yday[0][tm->month] + tm->day - 1);

    return ((wday % 7) + 7) % 7;
}

/* Compute the day of the year.  */
int dayOfYear(uritm *tm)
{
    return __mon_yday[IS_LEAP(tm->year)][tm->month] + (tm->day - 1);
}

#define ERRIF(x) \
    if (x) return 1;

int strptime_uri(const char *rp, uritm *tm)
{
    int lb, ub, nchar, val;

    ERRIF(!tm || !rp);

    /* Match year including century number. */
    lb    = 0;
    ub    = 9999;
    nchar = 4;
    val   = 0;

    while (*rp == '-' || *rp == ' ') ++rp;
    ERRIF(*rp < '0' || *rp > '9');

    do
    {
        val *= 10;
        val += *rp++ - '0';
    } while (--nchar > 0 && val * 10 <= ub && *rp >= '0' && *rp <= '9');

    ERRIF(val < lb || val > ub);

    tm->year = val;

    /* Match number of month. */
    ERRIF(!rp);
    lb    = 1;
    ub    = 12;
    nchar = 2;
    val   = 0;

    while (*rp == '-' || *rp == ' ') ++rp;
    ERRIF(*rp < '0' || *rp > '9');

    do
    {
        val *= 10;
        val += *rp++ - '0';
    } while (--nchar > 0 && val * 10 <= ub && *rp >= '0' && *rp <= '9');

    ERRIF(val < lb || val > ub);

    tm->month = val - 1;

    /* Match day of month. */
    ERRIF(!rp);
    lb    = 1;
    ub    = 31;
    nchar = 2;
    val   = 0;

    while (*rp == '-' || *rp == ' ') ++rp;
    ERRIF(*rp < '0' || *rp > '9');

    do
    {
        val *= 10;
        val += *rp++ - '0';
    } while (--nchar > 0 && val * 10 <= ub && *rp >= '0' && *rp <= '9');

    ERRIF(val < lb || val > ub);

    tm->day = val;

    return 0;
}

#undef ERRIF

int uritm_lessThan(uritm *tm1, uritm *tm2)
{
    if (!tm1 || !tm2) error("null pointer");
    if (tm1->year < tm2->year) return 1;
    if (tm1->month < tm2->month) return 1;
    if (tm1->day < tm2->day) return 1;
    return 0;
}

int uritm_greaterThan(uritm *tm1, uritm *tm2)
{
    if (!tm1 || !tm2) error("null pointer");

    if (tm1->year > tm2->year) return 1;
    if (tm1->month > tm2->month) return 1;
    if (tm1->day > tm2->day) return 1;
    return 0;
}

int uritm_equals(uritm *tm1, uritm *tm2)
{
    if (!tm1 || !tm2) error("Null pointer");

    return tm1->year == tm2->year && tm1->month == tm2->month &&
           tm1->day == tm2->day;
}

long uritm_difftime(uritm *tm1, uritm *tm2)
{
    long ans, i, c;
    if (!tm1 || !tm2) error("null pointer");

    ans  = dayOfYear(tm1) - dayOfYear(tm2);
    ans += 365 * (tm1->year - tm2->year);

    if (tm2->year < tm1->year)
    {
        for (i = tm2->year; i < tm1->year; ++i)
            if (IS_LEAP(i)) ++ans;

    } else if (tm1->year < tm2->year)
    {
        for (i = tm1->year; i < tm2->year; ++i)
            if (IS_LEAP(i)) --ans;
    }

    return ans;
}
